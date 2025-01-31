/*******************************************************************************
 * Copyright (c) 2023, 2024 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.vertiflex;

import java.util.List;

import org.eclipse.elk.alg.common.NodeMicroLayout;
import org.eclipse.elk.alg.vertiflex.options.VertiFlexOptions;
import org.eclipse.elk.alg.vertiflex.p1yplacement.NodeYPlacerStrategy;
import org.eclipse.elk.alg.vertiflex.p2relative.RelativeXPlacerStrategy;
import org.eclipse.elk.alg.vertiflex.p3absolute.AbsoluteXPlacerStrategy;
import org.eclipse.elk.alg.vertiflex.p4edgerouting.EdgerouterStrategy;
import org.eclipse.elk.core.AbstractLayoutProvider;
import org.eclipse.elk.core.UnsupportedConfigurationException;
import org.eclipse.elk.core.alg.AlgorithmAssembler;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.math.ElkMargin;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkNode;

/**
 * Layout provider for the y constraint tree layout algorithms.
 */
public final class VertiFlexLayoutProvider extends AbstractLayoutProvider {
    
    
    private final AlgorithmAssembler<VertiFlexLayoutPhases, ElkNode> algorithmAssembler =
        AlgorithmAssembler.<VertiFlexLayoutPhases, ElkNode>create(VertiFlexLayoutPhases.class);
    
    private double nodeNodeSpacing;

    @Override
    public void layout(final ElkNode graph, final IElkProgressMonitor progressMonitor) {
        List<ILayoutProcessor<ElkNode>> algorithm = assembleAlgorithm(graph);

        progressMonitor.begin("Tree layout", algorithm.size());
        
        nodeNodeSpacing = graph.getProperty(CoreOptions.SPACING_NODE_NODE);
        
        // if requested, compute nodes's dimensions, place node labels, ports, port labels, etc.
        if (!graph.getProperty(VertiFlexOptions.OMIT_NODE_MICRO_LAYOUT)) {
            NodeMicroLayout.forGraph(graph)
                           .execute();
        }
        
        // pre calculate the root node and save it
        ElkNode root = VertiFlexUtil.findRoot(graph);
        graph.setProperty(InternalProperties.ROOT_NODE, root);
        if (root == null) {
            throw new UnsupportedConfigurationException("The given graph is not a tree!");
        }
        
        for (ElkNode child : graph.getChildren()) {
            int numberOfParents;
            numberOfParents = child.getIncomingEdges().size();
            if (numberOfParents > 1) {
                throw new UnsupportedConfigurationException("The given graph is not an acyclic tree!");
            }
            
            // reset position
            child.setLocation(0, 0);
        }
        
        // check that vertical constraints are ordered in valid manner i.e. children always have higher vertical 
        // constraints than their parents
        checkVerticalConstraintValidity(root, 0);
        
        // store model order
        int count = 0;
        for (ElkNode node : graph.getChildren()) {
            node.setProperty(InternalProperties.NODE_MODEL_ORDER, count);
            count += 1;
        }

        for (ILayoutProcessor<ElkNode> processor : algorithm) {
            processor.process(graph, progressMonitor.subTask(1));
        }
        
        setGraphSize(graph);

        progressMonitor.done();
    }
    
    /**
     * Configure the layout provider by assembling different layout processors.
     * 
     * @param graph The graph which shall be laid out.
     * @return The list of assembled layout processors.
     */
    public List<ILayoutProcessor<ElkNode>> assembleAlgorithm(final ElkNode graph) {
        algorithmAssembler.reset();

        // Configure phases
        algorithmAssembler.setPhase(VertiFlexLayoutPhases.P1_NODE_Y_PLACEMENT,
                NodeYPlacerStrategy.SIMPLE_Y_PLACING);
        algorithmAssembler.setPhase(VertiFlexLayoutPhases.P2_NODE_RELATIVE_PLACEMENT,
                RelativeXPlacerStrategy.SIMPLE_X_PLACING);
        algorithmAssembler.setPhase(VertiFlexLayoutPhases.P3_NODE_ABSOLUTE_PLACEMENT,
                AbsoluteXPlacerStrategy.ABSOLUTE_XPLACING);
        
        EdgerouterStrategy routerStrategy;
        switch (graph.getProperty(VertiFlexOptions.LAYOUT_STRATEGY)) {
            case BEND:
                routerStrategy = EdgerouterStrategy.BEND_ROUTING;
                break;
            case STRAIGHT:
            default:
                routerStrategy = EdgerouterStrategy.DIRECT_ROUTING;
                break;
            
        }
        algorithmAssembler.setPhase(VertiFlexLayoutPhases.P4_EDGE_ROUTING, routerStrategy);

        // Assemble the algorithm
        return algorithmAssembler.build(graph);
    }
    
    /** Checks whether a vertical constraint is larger than the constraints set by any ancestor nodes.*/
    private void checkVerticalConstraintValidity(final ElkNode root, final double currentMinConstraint) {

        double rootHeight;
        if (root.hasProperty(VertiFlexOptions.VERTICAL_CONSTRAINT)) {
            rootHeight = root.getProperty(VertiFlexOptions.VERTICAL_CONSTRAINT);
        } else {
            rootHeight = currentMinConstraint;
        }
        
        double newMinConstraint = rootHeight + root.getHeight() 
            + Math.max(root.getProperty(CoreOptions.MARGINS).bottom, nodeNodeSpacing);
        
        for (ElkEdge outgoingEdge : root.getOutgoingEdges()) {
            ElkNode child = (ElkNode) outgoingEdge.getTargets().get(0);
            if (child.hasProperty(VertiFlexOptions.VERTICAL_CONSTRAINT)) {
                if (newMinConstraint > child.getProperty(VertiFlexOptions.VERTICAL_CONSTRAINT) 
                            + child.getProperty(CoreOptions.MARGINS).top) {
                    throw new UnsupportedConfigurationException("Invalid vertical constraints. Node "
                            + child.getIdentifier() + " has a vertical constraint that is too low for its ancestors.");
                }
            }
        }
        for (ElkEdge outgoingEdge : root.getOutgoingEdges()) {
            ElkNode child = (ElkNode) outgoingEdge.getTargets().get(0);
            checkVerticalConstraintValidity(child, newMinConstraint);
        }
    }
    
    /** Computes the space occupied by the layout and sets the graph size accordingly. */
    private void setGraphSize(final ElkNode graph) {
        ElkPadding padding = graph.getProperty(CoreOptions.PADDING);

        double maxX = 0.0;
        double maxY = 0.0;
        for (ElkNode node : graph.getChildren()) {
            ElkMargin margin = node.getProperty(CoreOptions.MARGINS);

            if (maxX < node.getX() + node.getWidth() + margin.right) {
                maxX = node.getX() + node.getWidth() + margin.right;
            }
            if (maxY < node.getY() + node.getHeight() + margin.bottom) {
                maxY = node.getY() + node.getHeight() + margin.bottom;
            }
        }

        graph.setWidth(maxX + padding.right);
        graph.setHeight(maxY + padding.bottom);
    }

}
