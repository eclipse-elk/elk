/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.

 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.disco;

import java.util.List;
import java.util.Set;

import org.eclipse.elk.alg.common.compaction.options.PolyominoOptions;
import org.eclipse.elk.alg.common.polyomino.PolyominoCompactor;
import org.eclipse.elk.alg.common.polyomino.structures.PlanarGrid;
import org.eclipse.elk.alg.common.polyomino.structures.Polyomino;
import org.eclipse.elk.alg.common.polyomino.structures.Polyominoes;
import org.eclipse.elk.alg.disco.graph.DCComponent;
import org.eclipse.elk.alg.disco.graph.DCGraph;
import org.eclipse.elk.alg.disco.options.DisCoOptions;
import org.eclipse.elk.alg.disco.structures.DCPolyomino;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.Quadruple;

import com.google.common.collect.Lists;

/**
 * <p>
 * Implementation of a connected components layout algorithm based on using polyominoes (i. e. low resolution
 * representations) derived from each connected component.
 * </p>
 * 
 * <p>
 * The approach is based on:
 * </p>
 * 
 * <p>
 * Freivalds, K., Dogrusoz, U. and Kikusts, P. Disconnected Graph Layout and the Polyomino Packing Approach. in Mutzel,
 * P., Jünger, M. and Leipert, S. (Eds.). Graph Drawing: 9th International Symposium, GD 2001 Vienna, Austria, September
 * 23-26, 2001 Revised Papers. Springer Berlin Heidelberg, Berlin, Heidelberg, 2002, 378-39. DOI: <a
 * href=http://dx.doi.org/10.1007/3-540-45848-4_30>10.1007/3-540-45848-4_30</a>
 * </p>
 * 
 * <p>
 * A subgroup of polyominoes might remind you of your sheltered childhood, when you still had time to play video games
 * on your popular, barely pocketable, japanese mobile video gaming device (see <a
 * href=https://en.wikipedia.org/wiki/Polyomino>Polyomino (Wikipedia)</a> and <a
 * href=https://en.wikipedia.org/wiki/Tetromino>Tetromino (Wikipedia)</a> for visual reminders, retrieved on
 * 2016-09-20).
 * </p>
 */
public class DisCoPolyominoCompactor implements ICompactor {

    ///////////////////////////////////////////////////////////////////////////////
    // Variables

    /**
     * number of spatial units of the {@link DCGraph} covered by the witdh (and height) of a discrete grid cell used for
     * creating {@link Polyomino polyominoes}.
     */
    private double gridCellSizeX, gridCellSizeY;
    /** upper bound on the size of a grid cell, 100 is taken as a recommendation from the paper. */
    private final double upperBound = 100.0;
    /** Graph to compact. */
    private DCGraph cmpGraph;
    /** {@link Polyomino Polyominoes} representing each one {@link DCComponent} of the given {@link DCGraph}. */
    private List<DCPolyomino> polys;
    /**
     * finite grid representing the infinite planar grid ℕ x ℕ on which the {@link DCPolyomino polyominoes} are placed.
     */
    private PlanarGrid grid;
    /**
     * Desired aspect ratio of the final layout (width by height).
     */
    private Double aspectRatio;

    /**
     * Determines whether holes in polyominoes are to be filled.
     */
    private boolean fill;

    ///////////////////////////////////////////////////////////////////////////////
    // Implemenation of interface method

    /*
     * (non-Javadoc)
     * 
     * @see org.eclipse.elk.alg.disco.ICompactor#compact(org.eclipse.elk.alg.disco.graph.DCGraph)
     */
    @Override
    public void compact(final DCGraph graph) {
        // Initialization
        this.cmpGraph = graph;
        polys = Lists.newArrayList();

        // Computation of the ideal step size of the discrete grid in relation to the original graphs resolution.
        double gridCellRecommendation = computeCellSize(cmpGraph);
        gridCellSizeX = gridCellRecommendation;
        gridCellSizeY = gridCellRecommendation;

        fill = cmpGraph.getProperty(PolyominoOptions.POLYOMINO_FILL);
        aspectRatio = cmpGraph.getProperty(CoreOptions.ASPECT_RATIO);
        if (aspectRatio == null) {
            aspectRatio = Double.valueOf(1.0);
        }

        if (aspectRatio > 1.0) {
            // Width greater than height
            gridCellSizeX *= aspectRatio;
        } else {
            gridCellSizeY /= aspectRatio;
        }

        // 1.) Convert DCComponents into polyominoes and generate a grid large enough to contain them all after
        // packing.
        createPolyominoes();
        // 2.) Do the actual layout of the polyominoes.
        packPolyominoes();

        // 3.) Apply layout back to the DCGraph
        applyToDCGraph();

        // Used for visual debugging purposes only.
        cmpGraph.setProperty(DisCoOptions.DEBUG_DISCO_POLYS, polys);
    }

    ///////////////////////////////////////////////////////////////////////////////
    // Private submethods

    /**
     * <p>
     * Computes the optimal step size for the polyomino grid, satisfying this equation from the paper.
     * </p>
     * 
     * <p>
     * (c*n-1)*l² - (Σ from i=1 to n: (W_i + H_i) * l) - (Σ from i=1 to n: W_i * H_i) = 0
     * </p>
     * 
     * <p>
     * with
     * <ul>
     * <li>c: upper bound on the avarage size of a polyomino (given as 100 in the paper based on experiments regarding
     * computation time and quality resulting layout)</li>
     * <li>l: grid step (result of this method)</li>
     * <li>n: number of connected components of the {@link DCGraph}</li>
     * <li>W_i: width of {@link DCComponent} i</li>
     * <li>H_i: height of {@link DCComponent} i</li>
     * </ul>
     * </p>
     * 
     * @param graph
     *            input graph to compute optimal step size for
     * @return optimal step size (number of spatial units of the {@link DCGraph} covered by the witdh or height of a
     *         discrete grid cell used for creating {@link DCPolyomino polyominoes})
     */
    private double computeCellSize(final DCGraph graph) {
        double sumTerm = 0.0;
        double prodTerm = 0.0;
        Set<DCComponent> comps = graph.getComponents();
        int numOfComps = comps.size();

        for (DCComponent comp : comps) {
            KVector bounds = comp.getDimensionsOfBoundingRectangle();
            double width = bounds.x;
            double height = bounds.y;
            sumTerm += width + height;
            prodTerm += width * height;
        }

        // Numerator of positive solution of quadratic equation (the negative solution is is not interesting for us, as
        // step size should be positive)
        final double four = 4.0;
        double numerator =
                Math.sqrt(four * upperBound * numOfComps * prodTerm - four * prodTerm + sumTerm * sumTerm) + sumTerm;
        double denominator = 2 * (upperBound * numOfComps - 1.0);

        // Shouldn't happen with reasonable values for parameter "bound", ignore denominator just in case
        if (denominator == 0.0) {
            return numerator;
        }
        return numerator / denominator;
    }

    /**
     * Creates {@link DCPolyomino polyominoes} from all {@link DCComponent DCComponents} of the {@link DCGraph} to
     * compact. Initializes an underlying grid structure to represent the polyomino based layout. Prerequisite:
     * {@link #computeCellSize(DCGraph)} must have been called before this method.
     */
    private void createPolyominoes() {
        DCPolyomino poly;

        Set<DCComponent> comps = cmpGraph.getComponents();

        for (DCComponent comp : comps) {
            poly = new DCPolyomino(comp, gridCellSizeX, gridCellSizeY);
            polys.add(poly);
        }
    }

    /**
     * Places {@link DCPolyomino polyominoes} close together, in order to achieve a minimum area for the final drawing,
     * based on the heuristic algorithm PackPolyominoes from the paper. Prerequisite: {@link #createPolyominoes()} must
     * have been called before this method.
     */
    private void packPolyominoes() {

        // Id is for debugging purposes only (s. DisCoDebugView).
        int id = 0;
        for (DCPolyomino poly : polys) {
            poly.setId(id);
            id++;
        }

        // Use the more generic polyomino compactor which isn't restricted to components
        Polyominoes<DCPolyomino> polyHolder = new Polyominoes<DCPolyomino>(polys, aspectRatio, fill);
        new PolyominoCompactor().packPolyominoes(polyHolder);

        polys = polyHolder.getPolyominoes();
        grid = polyHolder.getGrid();
    }

    /**
     * Sets an offset to the {@link DCComponent DCComponents} of the given graph indicating a vector, which yields its
     * new absolute position after packing, if added to its original coordinates. Prerequisite:
     * {@link #packPolyominoes()} must have been called before this method.
     */
    private void applyToDCGraph() {
        // Crop base grid to the bounding box of all filled/marked cells.
        Quadruple<Integer, Integer, Integer, Integer> gridCrop = grid.getFilledBounds();
        // Compute width and height of the minimum bounding box of the new layout and set the the dimensions of the
        // graph accordingly. The coordinates are thereby translated from the discrete integer coordinates of the low
        // resolution polyominoes back to the real coordinates of the DCGraph.
        ElkPadding padding = cmpGraph.getProperty(DisCoOptions.PADDING);
        double paddingHori = padding.getHorizontal();
        double paddingVert = padding.getVertical();
        double parentWidth = (gridCrop.getThird() * gridCellSizeX) + paddingHori;
        double parentHeight = (gridCrop.getFourth() * gridCellSizeY) + paddingVert;
        cmpGraph.setDimensions(new KVector(parentWidth, parentHeight));

        for (DCPolyomino poly : polys) {
            // Compute the position of the upper left corner of each polyomino relative to the upper left corner of the
            // bounding box of all filled/marked cells of the base grid.
            int absoluteIntPositionX = poly.getX() - gridCrop.getFirst();
            int absoluteIntPositionY = poly.getY() - gridCrop.getSecond();

            // Compute new position of DCComponent
            KVector absolutePositionOnCanvas = new KVector(absoluteIntPositionX, absoluteIntPositionY)
                    .scale(poly.getCellSizeX(), poly.getCellSizeY()).add(poly.getOffset());
            // Get original coordinates of DCCopmponent (These will not be updated, as so far all classes implementing
            // the IGraphTransformer only use the offsets of the components to layout their underlying graph structures.
            KVector originalCoordinates = poly.getRepresentee().getMinCorner();
            // Set offset of new position relative to the original position of the component
            poly.getRepresentee().setOffset(absolutePositionOnCanvas.sub(originalCoordinates));
        }
    }

}
