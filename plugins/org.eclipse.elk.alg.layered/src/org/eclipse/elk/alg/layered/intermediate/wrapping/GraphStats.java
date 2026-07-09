/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
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
package org.eclipse.elk.alg.layered.intermediate.wrapping;

import java.util.Iterator;
import java.util.List;
import java.util.function.BinaryOperator;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.options.Direction;

/**
 * Collects (lazily) some information about a layered graph. For instance: longest path and maximum width and height of
 * a layer.
 */
public class GraphStats {

    // CHECKSTYLEOFF Javadoc
    // CHECKSTYLEOFF VisibilityModifier
    public final LGraph graph;
    private double spacing = 0;
    private double inLayerSpacing = 0;

    // computed during instantiation
    public final double dar;
    public final int longestPath;
    
    // computed on demand 
    private Double maxWidth;
    private Double maxHeight;
    private Double sumWidth;
    
    private double[] widths;
    private double[] heights;
    
    private boolean[] cutsAllowed;

    /**
     * Initialize new stats object, most of the internal data are not yet initialized.
     */
    public GraphStats(final LGraph graph) {
        this.graph = graph;

        // since the graph direction may be horizontal (default) or vertical,
        //  the desired aspect ration must be adjusted correspondingly
        //  (ELK internally always assumes left to right, however the
        //  aspect ratio is not adjusted during graph import)
        final Direction dir = graph.getProperty(LayeredOptions.DIRECTION);
        final double aspectRatio = graph.getProperty(LayeredOptions.ASPECT_RATIO);
        final double correction = graph.getProperty(LayeredOptions.WRAPPING_CORRECTION_FACTOR);
        if (dir == Direction.LEFT || dir == Direction.RIGHT || dir == Direction.UNDEFINED) {
            this.dar = aspectRatio * correction;
        } else {
            this.dar = 1 / (aspectRatio * correction);
        }
        
        this.spacing = graph.getProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS);
        this.inLayerSpacing = graph.getProperty(LayeredOptions.SPACING_NODE_NODE);
        
        this.longestPath = graph.getLayers().size();
    }
    
    /* ------------------------------------------------------------------------------------------- */
    /* Approximations
    /* ------------------------------------------------------------------------------------------- */
   
    /** @return the max width of any node of the layer. */
    public double getApproximateLayerWidth(final Layer l) {
        return determineLayerWidth(l);
    }
    
    /** @return the approximate width of the <em>current</em> layering (the sum of each layer's width). */
    public double getApproximateLayeringWidth() {
        return getSumWidth();
    }

    /**
     * @return the approximate width of the <em>new</em> layering that would be created if the {@code cuts} were
     *         applied.
     */
    public double getApproximateChunkBasedLayeringWidth(final List<Integer> cuts) {
        if (cuts.isEmpty()) {
            return getApproximateLayeringWidth();
        }
        double width = 0;
        int lowIdx = 0;
        for (Integer highIdx : cuts) {
            width = Math.max(width, determineChunkWidth(lowIdx, highIdx));
            lowIdx = highIdx;
        }
        width = Math.max(width, determineChunkWidth(lowIdx, graph.getLayers().size()));
        return width;
    }
        
    /** @return the sum of the layer's node heights. */
    public double getApproximateLayerHeight(final Layer l) {
        return determineLayerHeight(l);
    }
    
    /** @return the approximate height of the <em>current</em> layering (the max height of any layer). */
    public double getApproximateLayeringHeight() {
        return getMaxHeight();
    }
    
    /**
     * @return the approximate height of the <em>new</em> layering that would be created if the {@code cuts} were
     *         applied.
     */
    public double getApproximateChunkBasedLayeringHeight(final List<Integer> cuts) {
        if (cuts.isEmpty()) {
            return 0;
        }
        double height = 0;
        int lowIdx = 0;
        for (Integer highIdx : cuts) {
            height += determineChunkHeight(lowIdx, highIdx);
            lowIdx = highIdx;
        }
        height += determineChunkHeight(lowIdx, graph.getLayers().size());
        return height;
    }
    
    /* ------------------------------------------------------------------------------------------- */
    /* Widths
    /* ------------------------------------------------------------------------------------------- */
        
    /**
     * @return the maxWidth
     */
    public double getMaxWidth() {
        if (maxWidth == null) {
            maxWidth = determineWidth(Math::max);
        }
        return maxWidth;
    }
    
    public double getSumWidth() {
        if (sumWidth == null) {
            sumWidth = determineWidth((a, b) -> a + b);
        }
        return sumWidth;
    }
    
    /**
     * @return the widths
     */
    public double[] getWidths() {
        if (widths == null) {
            initWidthsAndHeights();
        }
        return widths;
    }
    
    private void initWidthsAndHeights() {
        int n = longestPath;
        this.widths = new double[n];
        this.heights = new double[n];
        
        for (int i = 0; i < n; i++) {
            Layer l = graph.getLayers().get(i);
            widths[i] = determineLayerWidth(l);
            heights[i] = determineLayerHeight(l);
        }
    }
    
    private double determineWidth(final BinaryOperator<Double> fun) {
        return graph.getLayers().stream()
            .map(l -> determineLayerWidth(l))
            .reduce(fun).get();
    }
    
    /**
     * @return max of any node width in the layer plus spacing
     */
    private double determineLayerWidth(final Layer l) {
        double maxW = 0;
        for (LNode n : l.getNodes()) {
            double nW = n.getSize().x + n.getMargin().right + n.getMargin().left + spacing;
            maxW = Math.max(maxW, nW);
        }
        return maxW;
    }
    
    /**
     * @return the sum of layer widths between lowIdx (inclusive) and highIdx (exclusive).
     */
    private double determineChunkWidth(final int lowIdx, final int highIdx) {
        return graph.getLayers().subList(lowIdx, highIdx).stream()
                .mapToDouble(l -> determineLayerWidth(l))
                .reduce(Double::sum)
                .orElse(0);
    }

    /* ------------------------------------------------------------------------------------------- */
    /* Heights
    /* ------------------------------------------------------------------------------------------- */
    
    /**
     * @return the maxHeight
     */
    public double getMaxHeight() {
        if (maxHeight == null) {
            maxHeight = determineHeight(Math::max);
        }
        return maxHeight;
    }
    
    /**
     * @return the heights
     */
    public double[] getHeights() {
        if (heights == null) {
            initWidthsAndHeights();
        }
        return heights;
    }
    
    private double determineHeight(final BinaryOperator<Double> fun) {
        return graph.getLayers().stream()
            .map(l -> determineLayerHeight(l))
            .reduce(fun).get();
    }
    
    private double determineLayerHeight(final Layer layer) {
        double lH = 0;
        for (LNode n : layer.getNodes()) {
            
            lH += n.getSize().y + n.getMargin().bottom + n.getMargin().top + inLayerSpacing;
            
            for (LEdge inc : n.getIncomingEdges()) {
                if (inc.getSource().getNode().getType() == NodeType.NORTH_SOUTH_PORT) {
                    LNode src = inc.getSource().getNode();
                    LNode origin = (LNode) src.getProperty(InternalProperties.ORIGIN);
                    lH += origin.getSize().y + origin.getMargin().bottom + origin.getMargin().top; 
                }
            }
            
        }
        
        return lH;
    }
  
    /**
     * @return the maximum layer height between lowIdx (inclusive) and highIdx (exclusive).
     */
    private double determineChunkHeight(final int lowIdx, final int highIdx) {
        return graph.getLayers().subList(lowIdx, highIdx).stream()
                .mapToDouble(l -> determineLayerHeight(l))
                .max()
                .orElse(0);
    }
    
    /* ------------------------------------------------------------------------------------------- */
    /* Cutting allowed?
    /* ------------------------------------------------------------------------------------------- */

    /**
     * @param layerIndex
     *            the index of a layer.
     * @return whether it is allowed to but <em>before</em> {@code layerIndex}.
     */
    public boolean isCutAllowed(final int layerIndex) {
        if (cutsAllowed == null) {
            initCutAllowed();
        }
        return cutsAllowed[layerIndex];
    }
    
    private void initCutAllowed() {
        if (cutsAllowed != null) {
            return;
        }
        
        cutsAllowed = new boolean[graph.getLayers().size()];
        cutsAllowed[0] = false;
        
        if (graph.hasProperty(LayeredOptions.WRAPPING_VALIDIFY_FORBIDDEN_INDICES)) {
            // user-specified forbidden indices
            List<Integer> forbidden = graph.getProperty(LayeredOptions.WRAPPING_VALIDIFY_FORBIDDEN_INDICES);
            for (int f : forbidden) {
                if (f > 0 && f < cutsAllowed.length) {
                    cutsAllowed[f] = false;
                }
            }
        } else {
            // 'default' behavior as implemented by 'isCutAllowed(..)'
            Iterator<Layer> layerIt = graph.getLayers().iterator();
            // skip the first layer
            if (layerIt.hasNext()) {
                layerIt.next();
            }
            // now check for valid cut indexes
            int i = 1;
            while (layerIt.hasNext()) {
                Layer layer = layerIt.next();
                cutsAllowed[i++] = isCutAllowed(layer);
            }
        }
    }
    
    /**
     * @return the cutsAllowed
     */
    public boolean[] getCutsAllowed() {
        if (cutsAllowed == null) {
            initCutAllowed();
        }
        return cutsAllowed;
    }
    
    /**
     * @return whether it is allowed to cut <em>before</em> the given layer.
     */
    private boolean isCutAllowed(final Layer layer) {
        // we only allow to cut between a pair of layers 
        // if there is only one pair 
        // of nodes that is connected by 1 or more edges
        boolean cutAllowed = true;
        LNode n1 = null, n2 = null;
        check: 
        for (LNode tgt: layer.getNodes()) {
            for (LEdge e : tgt.getIncomingEdges()) {
                // check for different source
                if (n1 != null && n1 != tgt) {
                    cutAllowed = false;
                    break check;
                }
                n1 = tgt;
                // check for different target
                LNode src = e.getSource().getNode();
                if (n2 != null && n2 != src) {
                    cutAllowed = false;
                    break check;
                }
                n2 = src;
            }
        }
        return cutAllowed;
    }
    
}