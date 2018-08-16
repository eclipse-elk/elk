/*******************************************************************************
 * Copyright (c) 2013, 2015 Kiel University and others.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 *
 * Contributors:
 *     Kiel University - initial API and implementation
 *******************************************************************************/
package org.eclipse.elk.alg.mrtree.p4route;

import org.eclipse.elk.alg.mrtree.TreeLayoutPhases;
import org.eclipse.elk.alg.mrtree.graph.TEdge;
import org.eclipse.elk.alg.mrtree.graph.TGraph;
import org.eclipse.elk.alg.mrtree.graph.TNode;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * TODO: implement smart edge routing
 * 
 * This class implements a dull edge routing by setting just source and target of a edge.
 * 
 * @author sor
 * @author sgu
 * 
 */
public class EdgeRouter implements ILayoutPhase<TreeLayoutPhases, TGraph> {

    /** intermediate processing configuration. */
    private static final LayoutProcessorConfiguration<TreeLayoutPhases, TGraph> INTERMEDIATE_PROCESSING_CONFIG =
            LayoutProcessorConfiguration.<TreeLayoutPhases, TGraph>create();

    /**
     * {@inheritDoc}
     */
    @Override
    public LayoutProcessorConfiguration<TreeLayoutPhases, TGraph> getLayoutProcessorConfiguration(final TGraph graph) {
        return INTERMEDIATE_PROCESSING_CONFIG;
    }

    /**
     * {@inheritDoc}
     */
    public void process(final TGraph tGraph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Dull edge routing", 1);

        for (TNode tnode : tGraph.getNodes()) {
            for (TEdge tedge : tnode.getOutgoingEdges()) {
                tedge.getBendPoints().clear();
            }
        }
    }

}
