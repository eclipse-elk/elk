/*******************************************************************************
 * Copyright (c) 2023 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.vertiflex;

import java.util.List;

import org.eclipse.elk.alg.vertiflex.options.VertiFlexOptions;
import org.eclipse.elk.alg.vertiflex.p1yplacement.NodeYPlacerStrategy;
import org.eclipse.elk.alg.vertiflex.p2relative.RelativeXPlacerStrategy;
import org.eclipse.elk.alg.vertiflex.p3absolute.AbsoluteXPlacerStrategy;
import org.eclipse.elk.alg.vertiflex.p4edgerouting.EdgerouterStrategy;
import org.eclipse.elk.core.AbstractLayoutProvider;
import org.eclipse.elk.core.UnsupportedConfigurationException;
import org.eclipse.elk.core.alg.AlgorithmAssembler;
import org.eclipse.elk.core.alg.ILayoutProcessor;
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

    @Override
    public void layout(final ElkNode graph, final IElkProgressMonitor progressMonitor) {
        List<ILayoutProcessor<ElkNode>> algorithm = assembleAlgorithm(graph);

        progressMonitor.begin("Tree layout", algorithm.size());
        
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
        }
        
        // check that vertical constraints are ordered in valid manner i.e. children always have higher vertical 
        // constraints than their parents
        checkVerticalConstraintValidity(root);

        for (ILayoutProcessor<ElkNode> processor : algorithm) {
            processor.process(graph, progressMonitor.subTask(1));
        }

        progressMonitor.done();
    }
    
    /**
     * Configure the layout provider by assembling different layout processors.
     * 
     * @param graph The graph which shall be layout.
     * @return The list of assembled layout processors.
     */
    public List<ILayoutProcessor<ElkNode>> assembleAlgorithm(final ElkNode graph) {
        algorithmAssembler.reset();

        // Configure phases
        algorithmAssembler.setPhase(VertiFlexLayoutPhases.P1_NODE_Y_PLACEMENT,
                NodeYPlacerStrategy.SIMPLE_YPLACING);
        algorithmAssembler.setPhase(VertiFlexLayoutPhases.P2_NODE_RELATIVE_PLACEMENT,
                RelativeXPlacerStrategy.SIMPLE_XPLACING);
        algorithmAssembler.setPhase(VertiFlexLayoutPhases.P3_NODE_ABSOLUTE_PLACEMENT,
                AbsoluteXPlacerStrategy.ABSOLUTE_XPLACING);
        algorithmAssembler.setPhase(VertiFlexLayoutPhases.P4_EDGE_ROUTING,
                EdgerouterStrategy.DIRECT_ROUTING);

        // Assemble the algorithm
        return algorithmAssembler.build(graph);
    }
    
    private void checkVerticalConstraintValidity(final ElkNode root) {
        if (root.hasProperty(VertiFlexOptions.VERTICAL_CONSTRAINT)) {
            double rootHeight = root.getProperty(VertiFlexOptions.VERTICAL_CONSTRAINT);
            for (ElkEdge outgoingEdge : root.getOutgoingEdges()) {
                ElkNode child = (ElkNode) outgoingEdge.getTargets().get(0);
                if (child.hasProperty(VertiFlexOptions.VERTICAL_CONSTRAINT)) {
                    if (rootHeight + root.getHeight() + root.getProperty(CoreOptions.MARGINS).bottom
                            >= child.getProperty(VertiFlexOptions.VERTICAL_CONSTRAINT) 
                                + child.getProperty(CoreOptions.MARGINS).top) {
                        throw new UnsupportedConfigurationException("Invalid vertical constraints. Node " 
                            + root.getIdentifier() + " must have a smaller vertical constraint than its child " 
                                + child.getIdentifier() + ". This includes both node's margins.");
                    }
                }
            }
        }
        for (ElkEdge outgoingEdge : root.getOutgoingEdges()) {
            ElkNode child = (ElkNode) outgoingEdge.getTargets().get(0);
            checkVerticalConstraintValidity(child);
        }
    }

}
