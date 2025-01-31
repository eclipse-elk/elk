/*******************************************************************************
 * Copyright (c) 2022 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0  *******************************************************************************/
package org.eclipse.elk.alg.mrtree.intermediate;

import java.util.ArrayList;
import java.util.List;
import org.eclipse.elk.alg.mrtree.graph.TGraph;
import org.eclipse.elk.alg.mrtree.graph.TNode;
import org.eclipse.elk.alg.mrtree.options.InternalProperties;
import org.eclipse.elk.alg.mrtree.options.MrTreeOptions;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.core.util.Pair;

/**
 * A processor that computes the start and end coordinates for each levels nodes.
 * The used coordinates depend on the layout direction of the graph.
 * 
 * @author jnc, sdo
 */
public class LevelCoordinatesProcessor implements ILayoutProcessor<TGraph> {
    /**
     * {@inheritDoc}
     */
    public void process(final TGraph tGraph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Processor determine the coords for each level", 1f);

        List<Pair<Double, Double>> levels = new ArrayList<Pair<Double, Double>>();
        
        // Set up levels
        for (TNode n : tGraph.getNodes()) {
            while (n.getProperty(MrTreeOptions.TREE_LEVEL) > levels.size() - 1) {
                levels.add(new Pair<>(Double.MAX_VALUE, -Double.MAX_VALUE));
            }
            
            int curLevel = n.getProperty(MrTreeOptions.TREE_LEVEL);
            if (tGraph.getProperty(MrTreeOptions.DIRECTION).isHorizontal()) {
                if (n.getPosition().x < levels.get(curLevel).getFirst()) {
                    levels.get(curLevel).setFirst(n.getPosition().x);
                }
                if (n.getPosition().x + n.getSize().x > levels.get(curLevel).getSecond()) {
                    levels.get(curLevel).setSecond(n.getPosition().x + n.getSize().x);
                }
            } else {
                if (n.getPosition().y < levels.get(curLevel).getFirst()) {
                    levels.get(curLevel).setFirst(n.getPosition().y);
                }
                if (n.getPosition().y + n.getSize().y > levels.get(curLevel).getSecond()) {
                    levels.get(curLevel).setSecond(n.getPosition().y + n.getSize().y);
                }
            }
        }
        
        // Set node properties
        for (TNode n : tGraph.getNodes()) {
            int curLevel = n.getProperty(MrTreeOptions.TREE_LEVEL);
            n.setProperty(InternalProperties.LEVELMIN, levels.get(curLevel).getFirst());
            n.setProperty(InternalProperties.LEVELMAX, levels.get(curLevel).getSecond());
        }
        
        progressMonitor.done();
    }
}
