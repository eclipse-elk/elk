/*******************************************************************************
 * Copyright (c) 2022 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.mrtree.intermediate;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.eclipse.elk.alg.mrtree.graph.TGraph;
import org.eclipse.elk.alg.mrtree.graph.TNode;
import org.eclipse.elk.alg.mrtree.options.InternalProperties;
import org.eclipse.elk.alg.mrtree.options.MrTreeOptions;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * A processor that computes the treeLevel property for each node.
 * 
 * @author jnc, sdo
 */
public class LevelProcessor implements ILayoutProcessor<TGraph> {
    /** this map temporarily contains the level of each node in the given graph. */
    private Map<Integer, Integer> gloLevelMap = new HashMap<Integer, Integer>();
    
    /**
     * {@inheritDoc}
     */
    public void process(final TGraph tGraph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Processor compute fanout", 1);
        
        List<TNode> roots = tGraph.getNodes().stream().
                filter(x -> x.getProperty(InternalProperties.ROOT)).
                collect(Collectors.toList());
        
        setLevel(roots, 0);
        
        for (TNode tNode : tGraph.getNodes()) {
            int level = gloLevelMap.get(tNode.id) != null ? gloLevelMap.get(tNode.id) : 0;
            tNode.setProperty(MrTreeOptions.TREE_LEVEL, level);
        }

        progressMonitor.done();
    }

    /**
     * Sets the level for all nodes.
     *  
     * @param currentLevel the level we are currently in
     * @param level the level index
     */
    private void setLevel(final List<TNode> currentLevel, final int level) {
        if (!currentLevel.isEmpty()) {
            LinkedList<TNode> nextLevel = new LinkedList<TNode>();
            
            for (TNode tNode : currentLevel) {
                gloLevelMap.put(tNode.id, level);
                for (TNode tChild : tNode.getChildren()) {
                    nextLevel.add(tChild);
                }
            }
            
            setLevel(nextLevel, level + 1);
        }
    }
}
