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
package org.eclipse.elk.alg.spore.p1structure;

import java.util.List;

import org.eclipse.elk.alg.common.BowyerWatsonTriangulation;
import org.eclipse.elk.alg.common.spore.InternalProperties;
import org.eclipse.elk.alg.spore.SPOrEPhases;
import org.eclipse.elk.alg.spore.graph.Graph;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.core.util.IElkProgressMonitor;

import com.google.common.collect.Lists;

/**
 * This phase creates a Delaunay triangulation for the center vertices of the {@link Graph}.
 */
public class DelaunayTriangulationPhase implements ILayoutPhase<SPOrEPhases, Graph> {
    
    @Override
    public LayoutProcessorConfiguration<SPOrEPhases, Graph> getLayoutProcessorConfiguration(final Graph graph) {
        return LayoutProcessorConfiguration.<SPOrEPhases, Graph>create();
    }
    
    @Override
    public void process(final Graph graph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Delaunay triangulation", 1);
        
        List<KVector> vertices = Lists.newArrayList();
        graph.vertices.forEach(v -> vertices.add(v.originalVertex));
        
        String debugOutput = null;
        if (graph.getProperty(InternalProperties.DEBUG_SVG)) {
            debugOutput = ElkUtil.debugFolderPath("spore") + "10bw";
        }
        
        if (graph.tEdges == null) {
            graph.tEdges = BowyerWatsonTriangulation.triangulate(vertices, debugOutput);
        } else {
            graph.tEdges.addAll(BowyerWatsonTriangulation.triangulate(vertices, debugOutput));
        }

        progressMonitor.done();
    }
}
