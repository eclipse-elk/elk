/*******************************************************************************
 * Copyright (c) 2011, 2020 Kiel University and others.
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
package org.eclipse.elk.alg.force;

import java.util.List;
import java.util.Random;

import org.eclipse.elk.alg.common.NodeMicroLayout;
import org.eclipse.elk.alg.force.graph.FGraph;
import org.eclipse.elk.alg.force.model.AbstractForceModel;
import org.eclipse.elk.alg.force.model.EadesModel;
import org.eclipse.elk.alg.force.model.FruchtermanReingoldModel;
import org.eclipse.elk.alg.force.options.ForceModelStrategy;
import org.eclipse.elk.alg.force.options.ForceOptions;
import org.eclipse.elk.alg.force.options.InternalProperties;
import org.eclipse.elk.core.AbstractLayoutProvider;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkNode;

/**
 * Layout provider for the force layout algorithms.
 */
public final class ForceLayoutProvider extends AbstractLayoutProvider {

    /** the force model used by this layout algorithm. */
    private AbstractForceModel forceModel;
    /** connected components processor. */
    private ComponentsProcessor componentsProcessor = new ComponentsProcessor();

    @Override
    public void layout(final ElkNode elkGraph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("ELK Force", 1);
        
        // if requested, compute nodes's dimensions, place node labels, ports, port labels, etc.
        if (!elkGraph.getProperty(ForceOptions.OMIT_NODE_MICRO_LAYOUT)) {
            NodeMicroLayout.forGraph(elkGraph)
                           .execute();
        }
        
        // transform the input graph
        IGraphImporter<ElkNode> graphImporter = new ElkGraphImporter();
        FGraph fgraph = graphImporter.importGraph(elkGraph);

        // set special properties for the layered graph
        setOptions(fgraph, elkGraph);

        // update the force model depending on user selection
        updateModel(fgraph.getProperty(ForceOptions.MODEL));
        
        // split the input graph into components
        List<FGraph> components = componentsProcessor.split(fgraph);
        
        // perform the actual layout
        for (FGraph comp : components) {
            forceModel.layout(comp, progressMonitor.subTask(1.0f / components.size()));
        }
        
        // pack the components back into one graph
        fgraph = componentsProcessor.recombine(components);
        
        // apply the layout results to the original graph
        graphImporter.applyLayout(fgraph);
        
        progressMonitor.done();
    }
    
    /**
     * Set special layout options for the force graph.
     * 
     * @param fgraph a new force graph
     * @param parent the original parent node
     */
    private void setOptions(final FGraph fgraph, final ElkNode parent) {
        // set the random number generator based on the random seed option
        Integer randomSeed = fgraph.getProperty(ForceOptions.RANDOM_SEED);
        if (randomSeed != null) {
            int val = randomSeed;
            if (val == 0) {
                fgraph.setProperty(InternalProperties.RANDOM, new Random());
            } else {
                fgraph.setProperty(InternalProperties.RANDOM, new Random(val));
            }
        } else {
            fgraph.setProperty(InternalProperties.RANDOM, new Random(1));
        }
    }
    
    /**
     * Update the force model depending on the chosen strategy.
     * 
     * @param strategy a force model strategy
     */
    private void updateModel(final ForceModelStrategy strategy) {
        switch (strategy) {
        case EADES:
            if (!(forceModel instanceof EadesModel)) {
                forceModel = new EadesModel();
            }
            break;
        case FRUCHTERMAN_REINGOLD:
            if (!(forceModel instanceof FruchtermanReingoldModel)) {
                forceModel = new FruchtermanReingoldModel();
            }
            break;
        }
    }

}
