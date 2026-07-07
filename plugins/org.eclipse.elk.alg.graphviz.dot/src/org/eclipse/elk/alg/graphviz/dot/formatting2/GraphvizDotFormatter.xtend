/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.graphviz.dot.formatting2

import com.google.inject.Inject
import org.eclipse.elk.alg.graphviz.dot.dot.Attribute
import org.eclipse.elk.alg.graphviz.dot.dot.AttributeStatement
import org.eclipse.elk.alg.graphviz.dot.dot.EdgeStatement
import org.eclipse.elk.alg.graphviz.dot.dot.EdgeTarget
import org.eclipse.elk.alg.graphviz.dot.dot.Graph
import org.eclipse.elk.alg.graphviz.dot.dot.GraphvizModel
import org.eclipse.elk.alg.graphviz.dot.dot.Node
import org.eclipse.elk.alg.graphviz.dot.dot.NodeStatement
import org.eclipse.elk.alg.graphviz.dot.dot.Port
import org.eclipse.elk.alg.graphviz.dot.dot.Subgraph
import org.eclipse.elk.alg.graphviz.dot.services.GraphvizDotGrammarAccess
import org.eclipse.xtext.formatting2.AbstractFormatter2
import org.eclipse.xtext.formatting2.IFormattableDocument

import static org.eclipse.elk.alg.graphviz.dot.dot.DotPackage.Literals.*

class GraphvizDotFormatter extends AbstractFormatter2 {

    @Inject extension GraphvizDotGrammarAccess

    def dispatch void format(GraphvizModel graphvizModel, extension IFormattableDocument document) {
        var region = graphvizModel.previousHiddenRegion
        while (region !== null) {
            document.set(region, [autowrap])
            region = region.nextHiddenRegion
        }
        for (graph : graphvizModel.graphs) {
            graph.append[newLine]
            graph.format()
        }
    }

    def dispatch void format(Graph graph, extension IFormattableDocument document) {
        graph.regionFor.feature(GRAPH__TYPE).surround[oneSpace]
        interior(
            graph.regionFor.keyword('{').prepend[oneSpace].append[newLine],
            graph.regionFor.keyword('}'),
            [indent]
        )
        for (statement : graph.statements) {
            statement.append[newLine]
            statement.format()
        }
    }
    
    def dispatch void format(Attribute attribute, extension IFormattableDocument document) {
        attribute.regionFor.keyword(';').prepend[noSpace]
        attribute.regionFor.keyword('=').surround[oneSpace]
    }
    
    def dispatch void format(NodeStatement nodeStatement, extension IFormattableDocument document) {
        nodeStatement.regionFor.keyword(';').prepend[noSpace]
        nodeStatement.regionFor.keyword('[').prepend[oneSpace].append[noSpace]
        nodeStatement.regionFor.keyword(']').prepend[noSpace]
        nodeStatement.regionFor.keywords(',').forEach[prepend[noSpace].append[oneSpace]]
        nodeStatement.node?.format()
        for (attribute : nodeStatement.attributes) {
            attribute.format()
        }
    }
    
    def dispatch void format(Node node, extension IFormattableDocument document) {
        node.regionFor.keyword(':').surround[noSpace]
        node.port?.format()
    }
    
    def dispatch void format(EdgeStatement edgeStatement, extension IFormattableDocument document) {
        edgeStatement.regionFor.keyword(';').prepend[noSpace]
        edgeStatement.regionFor.keyword('[').prepend[oneSpace].append[noSpace]
        edgeStatement.regionFor.keyword(']').prepend[noSpace]
        edgeStatement.regionFor.keywords(',').forEach[prepend[noSpace].append[oneSpace]]
        edgeStatement.sourceNode?.format()
        for (target : edgeStatement.edgeTargets) {
            target.prepend[oneSpace]
            target.format()
        }
        for (attribute : edgeStatement.attributes) {
            attribute.format()
        }
    }
    
    def dispatch void format(EdgeTarget edgeTarget, extension IFormattableDocument document) {
        edgeTarget.regionFor.feature(EDGE_TARGET__OPERATOR).append[oneSpace]
        edgeTarget.targetSubgraph?.format()
        edgeTarget.targetnode?.format()
    }
    
    def dispatch void format(AttributeStatement attributeStatement, extension IFormattableDocument document) {
        attributeStatement.regionFor.keyword(';').prepend[noSpace]
        attributeStatement.regionFor.keyword('[').prepend[oneSpace].append[noSpace]
        attributeStatement.regionFor.keyword(']').prepend[noSpace]
        attributeStatement.regionFor.keywords(',').forEach[prepend[noSpace].append[oneSpace]]
        for (attribute : attributeStatement.attributes) {
            attribute.format()
        }
    }
    
    def dispatch void format(Subgraph subgraph, extension IFormattableDocument document) {
        subgraph.regionFor.keyword(';').prepend[noSpace]
        subgraph.regionFor.keyword(subgraphAccess.subgraphKeyword_1_0).append[oneSpace]
        interior(
            subgraph.regionFor.keyword('{').prepend[oneSpace].append[newLine],
            subgraph.regionFor.keyword('}'),
            [indent]
        )
        for (statement : subgraph.statements) {
            statement.append[newLine]
            statement.format()
        }
    }
    
    def dispatch void format(Port port, extension IFormattableDocument document) {
        port.regionFor.keyword(':').surround[noSpace]
    }
    
}
