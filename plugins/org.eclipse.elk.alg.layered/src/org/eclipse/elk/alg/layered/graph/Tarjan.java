/*******************************************************************************
 * Copyright (c) 2025 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 
 *******************************************************************************/
package org.eclipse.elk.alg.layered.graph;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.Stack;

import org.eclipse.elk.alg.layered.options.InternalProperties;

/**
 * Tarjan implementation to be used during layered layout.
 */


public class Tarjan {
    
    public Tarjan (List<LEdge> edgesToBeReversed, List<Set<LNode>> stronglyConnectedComponents,
            HashMap <LNode,Integer> nodeToSCCID) {
        this.edgesToBeReversed = edgesToBeReversed;
        this.stronglyConnectedComponents = stronglyConnectedComponents;
        this.nodeToSCCID = nodeToSCCID;
    }
    
    private List<LEdge> edgesToBeReversed;
    private int index = 0;
    protected List<Set<LNode>> stronglyConnectedComponents; // FIXME Why no ordered set here? this is bad
    private Stack<LNode> stack = new Stack<LNode>();
    private HashMap <LNode,Integer> nodeToSCCID = new HashMap<>();

    public void tarjan(final LGraph graph) {
        index = 0;
        stack = new Stack<LNode>();
        for (LNode node : graph.getLayerlessNodes()) {
            if (node.getProperty(InternalProperties.TARJAN_ID) == -1) {
                stronglyConnected(node);
                stack.clear();
            }
        }
    }

    public void stronglyConnected(final LNode v) {
        v.setProperty(InternalProperties.TARJAN_ID, index);
        v.setProperty(InternalProperties.TARJAN_LOWLINK, index);
        index++;
        stack.push(v);
        v.setProperty(InternalProperties.TARJAN_ON_STACK, true);
        for (LEdge edge : v.getConnectedEdges()) {
            if (edge.getSource().getNode() != v && !edgesToBeReversed.contains(edge)) {
                continue;
            }
            if (edge.getSource().getNode() == v && edgesToBeReversed.contains(edge)) {
                continue;
            }
            LNode target = null;
            if (edge.getTarget().getNode() == v) {
                target = edge.getSource().getNode();
            } else {
                target = edge.getTarget().getNode();
            }
            if (target.getProperty(InternalProperties.TARJAN_ID) == -1) {
                stronglyConnected(target);
                v.setProperty(InternalProperties.TARJAN_LOWLINK, 
                        Math.min(v.getProperty(InternalProperties.TARJAN_LOWLINK),
                        target.getProperty(InternalProperties.TARJAN_LOWLINK)));
            } else if (target.getProperty(InternalProperties.TARJAN_ON_STACK)) {
                v.setProperty(InternalProperties.TARJAN_LOWLINK, 
                        Math.min(v.getProperty(InternalProperties.TARJAN_LOWLINK),
                                target.getProperty(InternalProperties.TARJAN_ID)));
            }
        }
        if (v.getProperty(InternalProperties.TARJAN_LOWLINK) == v.getProperty(InternalProperties.TARJAN_ID)) {
            Set<LNode> sCC = new HashSet<LNode>();
            LNode n = null;
            do {
                n = stack.pop();
                n.setProperty(InternalProperties.TARJAN_ON_STACK, false);
                sCC.add(n);
            } while (v != n);
            if (sCC.size() >1) {
                int index = stronglyConnectedComponents.size();
                stronglyConnectedComponents.add(sCC);
                for (LNode node : sCC) {
                    nodeToSCCID.put(node, index);
                }
            }
        }
    }

    public void resetTarjan(final LGraph graph) {
        for (LNode n : graph.getLayerlessNodes()) {
            n.setProperty(InternalProperties.TARJAN_ON_STACK, false);
            n.setProperty(InternalProperties.TARJAN_LOWLINK, -1);
            n.setProperty(InternalProperties.TARJAN_ID, -1);
            stack.clear();
            for (LEdge e : n.getConnectedEdges()) {
                e.setProperty(InternalProperties.IS_PART_OF_CYCLE, false);
            }
        }
    }
}
