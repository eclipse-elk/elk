/*******************************************************************************
 * Copyright (c) 2022 Kiel University and others.
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
package org.eclipse.elk.alg.layered.p3order;

import org.eclipse.elk.alg.layered.LayeredPhases;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.intermediate.IntermediateProcessorStrategy;
import org.eclipse.elk.alg.layered.options.GreedySwitchType;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * This class does nothing instead of minimizing crossing. Everything remains as it was before.
 * <dl>
 *  <dt>Precondition:</dt>
 *      <dd>The graph has a proper layering, i.e. all long edges have been splitted; all nodes have at least fixed port
 *          sides. This is, however, not checked but preceding processors or phases might depend on this.</dd>
 *  <dt>Postcondition:</dt>
 *      <dd>The order of nodes and ports remains as it was before this phase.</dd>
 *  <dt>Dependencies:</dt>
 *      <dd>It is advised to set {@link GreedySwitchType} to OFF. Otherwise, nodes and ports might be reordered
 *      if crossings are created.</dd>
 * </dl>
 *
 */
public class NoCrossingMinimizer implements ILayoutPhase<LayeredPhases, LGraph> {

    @Override
    public void process(final LGraph graph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("No crossing minimization", 1);
        progressMonitor.done();
        
    }

    @Override
    public LayoutProcessorConfiguration<LayeredPhases, LGraph> getLayoutProcessorConfiguration(final LGraph graph) {
        LayoutProcessorConfiguration<LayeredPhases, LGraph> configuration =
                LayoutProcessorConfiguration.createFrom(INTERMEDIATE_PROCESSING_CONFIGURATION);
        configuration.addBefore(LayeredPhases.P3_NODE_ORDERING, IntermediateProcessorStrategy.PORT_LIST_SORTER);

        return configuration;
    }
    


    /** intermediate processing configuration. */
    private static final LayoutProcessorConfiguration<LayeredPhases, LGraph> INTERMEDIATE_PROCESSING_CONFIGURATION =
            LayoutProcessorConfiguration.<LayeredPhases, LGraph>create()
                    .addBefore(LayeredPhases.P3_NODE_ORDERING, IntermediateProcessorStrategy.LONG_EDGE_SPLITTER)
                    .addBefore(LayeredPhases.P4_NODE_PLACEMENT,
                            IntermediateProcessorStrategy.IN_LAYER_CONSTRAINT_PROCESSOR)
                    .after(LayeredPhases.P5_EDGE_ROUTING)
                        .add(IntermediateProcessorStrategy.LONG_EDGE_JOINER);

}
