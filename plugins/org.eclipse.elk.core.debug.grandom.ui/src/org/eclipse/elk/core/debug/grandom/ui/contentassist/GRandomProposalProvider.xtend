/*******************************************************************************
 * Copyright (c) 2018 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.debug.grandom.ui.contentassist

import com.google.common.collect.ImmutableMap
import org.eclipse.elk.core.debug.grandom.gRandom.Configuration
import org.eclipse.elk.core.debug.grandom.gRandom.Form
import org.eclipse.elk.core.debug.grandom.gRandom.RandGraph
import org.eclipse.emf.ecore.EObject
import org.eclipse.jface.viewers.StyledString
import org.eclipse.xtext.Keyword
import org.eclipse.xtext.ui.editor.contentassist.ContentAssistContext
import org.eclipse.xtext.ui.editor.contentassist.ICompletionProposalAcceptor

/**
 * See https://www.eclipse.org/Xtext/documentation/304_ide_concepts.html#content-assist
 * on how to customize the content assistant.
 */
class GRandomProposalProvider extends AbstractGRandomProposalProvider {

    static val notTogether = newArrayList(
        Pair.of([Keyword k | k.value == "maxWidth"], [Configuration c | c.form != Form.TREES]),
        Pair.of([Keyword k | k.value == "maxDegree"], [Configuration c | c.form != Form.TREES]),
        Pair.of([Keyword k | k.value == "partitionFraction"], [Configuration c | c.form != Form.BIPARTITE])
    )

    static val documentation = ImmutableMap.builder
        .put("nodes", "nodes".n.o.d(30))
        .put("edges", "edges".n.d("total"))
        .put("density", "Fraction of number of nodes squared.".o)
        .put("relative", "Relative to number of nodes.".o)
        .put("outgoing", "Number of outgoing edges per node.".o)
        .put("total", "Total number of edges in graph.".o.d(20))
        .put("labels", "Add labels.")
        .put("self loops", "Allow edges with equal source and target.")
        .put("filename", "Name + index.".d("\"random\""))
        .put("format", "Textual or xml.".d("kgt"))
        .put("hierarchy", "".e)
        .put("maxDegree", "Maximum Degree (trees)")
        .put("maxWidth", "Maximum Width (trees)")
        .put("partitionFraction", "Minimal fraction of nodes in second partition set (bipartite).")
        .put("seed", "Random seed.")
        .put("size", "size".o.d("node: width, height 30, port: width, height: 4"))
        .put("ports", "ports".o)
        .put("constraint", "set port constraints".o.d("free"))
        .put("incoming", "incoming ports".n.rel)
        // TODO needs separate default documentation.
        .put("north", "north ports".n.rel)
        .put("east", "east ports".n.rel)
        .put("south", "south ports".n.rel)
        .put("west", "west ports".n.rel)
        .put("re-use", "Fraction of edges with same source or target port.")
        .put("+/-", "Gaussian distribution: Mean +/- Standard deviation.")
        .put("to", "Equal distribution: Min to max value.")
        .build
    
    def static <T> String d(String s, T defaultVal){
        s + " (Default: " + defaultVal + ")"
    }
    def static String rel(String s){
        s + ". Relative to other given numbers."
    }
    def static String n(String s){
        "Number of " + s
    }
    def static String e(String s){
        s + "Experimental. Heare bee dragonns"
    }

    def static String o(String s){
        s + " With extra options in own block, i.e. {...}."
    }
    
    override StyledString getKeywordDisplayString(Keyword keyword) {
        val value = keyword.value
        val styled = new StyledString(value)
        val doc = documentation.get(value)
        if (doc !== null) 
            styled.append(new StyledString(": " + doc, StyledString.COUNTER_STYLER))
        return styled
    }
    
    override completeKeyword(Keyword keyword, ContentAssistContext contentAssistContext, 
        ICompletionProposalAcceptor acceptor){
        val conf = findConfig(contentAssistContext.rootModel as RandGraph, contentAssistContext.currentModel)
        for (rule : notTogether){
            if (rule.key.apply(keyword) && rule.value.apply(conf))
                return
        }
        super.completeKeyword(keyword, contentAssistContext, acceptor)
    }
    
    def Configuration findConfig(RandGraph root, EObject node) {
        if (root !== null && root.configs !== null)
            for (Configuration c: root.configs)
                if (node.in(c))
                    return c    
    }
    
    def boolean in(EObject key, EObject tree){
        for (v : tree.eContainer.eAllContents.toIterable) 
            if (key == v)
                return true
        return false
    }
    
}
