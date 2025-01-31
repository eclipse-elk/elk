/*******************************************************************************
 * Copyright (c) 2018, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.rectpacking;

import java.util.List;

import org.eclipse.elk.alg.common.NodeMicroLayout;
import org.eclipse.elk.alg.rectpacking.intermediate.IntermediateProcessorStrategy;
import org.eclipse.elk.alg.rectpacking.options.InternalProperties;
import org.eclipse.elk.alg.rectpacking.options.RectPackingOptions;
import org.eclipse.elk.core.AbstractLayoutProvider;
import org.eclipse.elk.core.alg.AlgorithmAssembler;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.eclipse.elk.core.util.BoxLayoutProvider;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkNode;

/**
 * A layout algorithm that does not take edges into account, but treats all nodes as isolated boxes. This is useful for
 * parts of a diagram that consist of objects without connections, such as parallel regions in Statecharts.
 * <p>
 * Nodes are viewed as rectangles and so {@link ElkNode}s are referred to as rectangles in the comments.
 * </p>
 * <p>
 * Depending on the settings, checks for a specified special case, calculates a layout with a approximation algorithm or
 * uses that approximation algorithm for the needed area of the rectangles and places the rectangles nicely aligned on
 * the drawing area according to that approximation.
 * </p>
 */
public class RectPackingLayoutProvider extends AbstractLayoutProvider {
    /** The thing that will assemble the processors that will constitute our algorithm. */
    private final AlgorithmAssembler<RectPackingLayoutPhases, ElkNode> algorithmAssembler =
            AlgorithmAssembler.<RectPackingLayoutPhases, ElkNode>create(RectPackingLayoutPhases.class);
    
    /**
     * Calculating and applying layout to the model.
     */
    @Override
    public void layout(final ElkNode layoutGraph, final IElkProgressMonitor progressMonitor) {
        
        progressMonitor.begin("Rectangle Packing", 1);

        ElkPadding padding = layoutGraph.getProperty(RectPackingOptions.PADDING);
        boolean fixedGraphSize = layoutGraph.getProperty(RectPackingOptions.NODE_SIZE_FIXED_GRAPH_SIZE);
        double nodeNodeSpacing = layoutGraph.getProperty(RectPackingOptions.SPACING_NODE_NODE);
        boolean tryBox = layoutGraph.getProperty(RectPackingOptions.TRYBOX);
        List<ElkNode> rectangles = layoutGraph.getChildren();
        
        // if requested, compute nodes's dimensions, place node labels, ports, port labels, etc.
        if (!layoutGraph.getProperty(RectPackingOptions.OMIT_NODE_MICRO_LAYOUT)) {
            NodeMicroLayout.forGraph(layoutGraph).execute();
        }
        
        // Check whether regions are stackable and do box layout instead.
        boolean stackable = false;
        if (tryBox && rectangles.size() >= 3) {
            ElkNode region1;
            ElkNode region2 = rectangles.get(0);
            ElkNode region3 = rectangles.get(1);
            int counter = 0;
            while (counter + 2 < rectangles.size()) {
                region1 = region2;
                region2 = region3;
                region3 = rectangles.get(counter + 2);
                if (region1.getHeight() >= region2.getHeight() + region3.getHeight() + nodeNodeSpacing
                        || region3.getHeight() >= region1.getHeight() + region2.getHeight() + nodeNodeSpacing) {
                    stackable = true;
                    break;
                } else {
                    counter++;
                }
            }
        } else {
            stackable = true;
        }
        if (!stackable) {
            // Set priority to invoke box layout.
            int priority = rectangles.size();
            for (ElkNode elkNode : rectangles) {
                elkNode.setProperty(CoreOptions.PRIORITY, priority);
                priority--;
            }
            new BoxLayoutProvider().layout(layoutGraph, new BasicProgressMonitor());
            progressMonitor.done();
            return;
        }

        List<ILayoutProcessor<ElkNode>> algorithm = assembleAlgorithm(layoutGraph);
        float monitorProgress = 1.0f / algorithm.size();
        
        // Invoke each layout processor
        int slotIndex = 0;
        for (ILayoutProcessor<ElkNode> processor : algorithm) {
            if (progressMonitor.isCanceled()) {
                return;
            }
            
            // Output debug graph
            // elkjs-exclude-start
            if (progressMonitor.isLoggingEnabled()) {
                progressMonitor.logGraph(layoutGraph, slotIndex + "-Before " + processor.getClass().getSimpleName());
            }
            // elkjs-exclude-end
            
            processor.process(layoutGraph, progressMonitor.subTask(monitorProgress));
            slotIndex++;
        }

        // Graph debug output
        // elkjs-exclude-start
        if (progressMonitor.isLoggingEnabled()) {
            progressMonitor.logGraph(layoutGraph, slotIndex + "-Finished");
        }
        // elkjs-exclude-end
        
        // Content alignment
        double realWidth = 0;
        double realHeight = 0;
        for (ElkNode rect : rectangles) {
            realWidth = Math.max(realWidth, rect.getX() + rect.getWidth());
            realHeight = Math.max(realHeight, rect.getY() + rect.getHeight());
        }

        ElkUtil.translate(layoutGraph,
                new KVector(
                        layoutGraph.getProperty(InternalProperties.DRAWING_WIDTH),
                        layoutGraph.getProperty(InternalProperties.DRAWING_HEIGHT)),
                new KVector(realWidth, realHeight));

        // Final touch.
        applyPadding(rectangles, padding);
        
        if (!fixedGraphSize) {
            ElkUtil.resizeNode(layoutGraph,
                    layoutGraph.getProperty(InternalProperties.DRAWING_WIDTH) + padding.getHorizontal(),
                    layoutGraph.getProperty(InternalProperties.DRAWING_HEIGHT) + padding.getVertical(), false, true);
        }

        // Do micro layout again since the whitspace elimination and other things might have changed node sizes.
        if (!layoutGraph.getProperty(RectPackingOptions.OMIT_NODE_MICRO_LAYOUT)) {
            NodeMicroLayout.forGraph(layoutGraph).execute();
        }
        // elkjs-exclude-start
        if (progressMonitor.isLoggingEnabled()) {
            progressMonitor.logGraph(layoutGraph, "Output");
        }
        // elkjs-exclude-end
        progressMonitor.done();
    }

    /**
     * Shifts all rectangles to the right and bottom according to the specified padding.
     * 
     * @param rectangles
     *            list of rectangles that have been placed.
     */
    private static void applyPadding(final List<ElkNode> rectangles, ElkPadding padding) {
        for (ElkNode rect : rectangles) {
            rect.setLocation(rect.getX() + padding.getLeft(), rect.getY() + padding.getTop());
        }
    }
    


    /**
     * Configure the layout provider by assembling different layout processors.
     * 
     * @param layoutGraph The graph which shall be layouted.
     * @return The list of assembled layout processors.
     */
    private List<ILayoutProcessor<ElkNode>> assembleAlgorithm(final ElkNode layoutGraph) {
        this.algorithmAssembler.reset();
        algorithmAssembler.setPhase(RectPackingLayoutPhases.P1_WIDTH_APPROXIMATION,
                layoutGraph.getProperty(RectPackingOptions.WIDTH_APPROXIMATION_STRATEGY));
        algorithmAssembler.setPhase(RectPackingLayoutPhases.P2_PACKING,
                layoutGraph.getProperty(RectPackingOptions.PACKING_STRATEGY));
        algorithmAssembler.setPhase(RectPackingLayoutPhases.P3_WHITESPACE_ELIMINATION,
                layoutGraph.getProperty(RectPackingOptions.WHITE_SPACE_ELIMINATION_STRATEGY));
        
        algorithmAssembler.addProcessorConfiguration(getPhaseIndependentLayoutProcessorConfiguration(layoutGraph));
        List<ILayoutProcessor<ElkNode>> processors = algorithmAssembler.build(layoutGraph);
        return processors;
    }
    

    
    /**
     * Returns an intermediate processing configuration with processors not tied to specific phases.
     * 
     * @param lgraph the layered graph to be processed. The configuration may vary depending on certain
     *               properties of the graph.
     * @return intermediate processing configuration. May be {@code null}.
     */
    private LayoutProcessorConfiguration<RectPackingLayoutPhases, ElkNode> getPhaseIndependentLayoutProcessorConfiguration(
            final ElkNode layoutGraph) {

        // Basic configuration
        LayoutProcessorConfiguration<RectPackingLayoutPhases, ElkNode> configuration =
                LayoutProcessorConfiguration.create();
        
        configuration.addBefore(RectPackingLayoutPhases.P1_WIDTH_APPROXIMATION,
                IntermediateProcessorStrategy.MIN_SIZE_PRE_PROCESSOR);
        
        configuration.addBefore(RectPackingLayoutPhases.P2_PACKING,
                IntermediateProcessorStrategy.MIN_SIZE_POST_PROCESSOR);
        
        if (layoutGraph.getProperty(RectPackingOptions.ORDER_BY_SIZE)) {
            configuration.addBefore(RectPackingLayoutPhases.P1_WIDTH_APPROXIMATION,
                    IntermediateProcessorStrategy.NODE_SIZE_REORDERER);
        }
        
        if (layoutGraph.getProperty(RectPackingOptions.INTERACTIVE)) {
            configuration.addBefore(RectPackingLayoutPhases.P1_WIDTH_APPROXIMATION,
                    IntermediateProcessorStrategy.INTERACTIVE_NODE_REORDERER);
        }
        
        return configuration;
    }
}
