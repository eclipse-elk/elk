/*******************************************************************************
 * Copyright (c) 2017 University of Kiel and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.common.polyomino.structures;

import java.util.List;

import org.eclipse.elk.alg.common.polyomino.ProfileFill;
import org.eclipse.elk.graph.properties.MapPropertyHolder;

import com.google.common.collect.Lists;

/**
 * Holds a bunch of polyominos and a grid to place them on.
 *
 * @param <P>
 *            a {@link Polyomino} or one of its derived classes
 */
public class Polyominoes<P extends Polyomino> extends MapPropertyHolder {
    private static final long serialVersionUID = -6003285848848775273L;
    private List<P> polys = Lists.newArrayList();
    private PlanarGrid grid;

    /**
     * Initializes the structure with some polyominoes and a desired aspect ratio for the underlying grid structure.
     * 
     * @param polys
     *            polyominoes or derivatives to add
     * @param aspectRatio
     *            desired aspect ratio given as the fraction width by height
     */
    public Polyominoes(final Iterable<P> polys, final double aspectRatio) {
        this(polys, aspectRatio, false);
    }

    /**
     * Initializes the structure with some polyominoes and a desired aspect ratio for the underlying grid structure.
     * Optionally holes in polyominoes can be filled usign{@link ProfileFill}.
     * 
     * @param polys
     *            polyominoes or derivatives to add
     * @param aspectRatio
     *            desired aspect ratio given as the fraction width by height
     * @param fill
     *            true, fill holes in polyominoes, false, keep polyominoes as they are
     */
    public Polyominoes(final Iterable<P> polys, final double aspectRatio, final boolean fill) {
        int gridWidth = 0;
        int gridHeight = 0;

        for (P poly : polys) {
            if (fill) {
                ProfileFill.fillPolyomino(poly);
            }
            this.polys.add(poly);
            gridWidth += poly.getWidth();
            gridHeight += poly.getHeight();
        }

        // Add width and height of the future center polyomino once again, for avoiding
        // Indices Out of bounds (might get better, see last todo)
        if (this.polys.size() > 0) {
            Polyomino poly = this.polys.get(0);
            gridWidth += poly.getWidth();
            gridHeight += poly.getHeight();
        }

        // Make the area 4 times greater for the case that all polyominoes can only be placed in one quadrant of the
        // coordinate system
        gridWidth *= 2;
        gridHeight *= 2;

        // Adjust grid size for desired aspect ratio
        if (aspectRatio > 1.0) {
            gridWidth = (int) Math.ceil(gridWidth * aspectRatio);
        } else {
            gridHeight = (int) Math.ceil(gridHeight / aspectRatio);
        }

        grid = new PlanarGrid(gridWidth, gridHeight);
    }

    /**
     * Returns the list of polyominoes.
     * 
     * @return list of polyominoes
     */
    public List<P> getPolyominoes() {
        return polys;
    }

    /**
     * Returns the underlying grid structure.
     * 
     * @return the grid structure the polyominoes can be placed on
     */
    public PlanarGrid getGrid() {
        return grid;
    }
}
