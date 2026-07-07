/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.graph.json.text.formatting2

import com.google.inject.Inject
import java.util.List
import org.eclipse.elk.graph.ElkEdge
import org.eclipse.elk.graph.ElkGraphElement
import org.eclipse.elk.graph.ElkLabel
import org.eclipse.elk.graph.ElkNode
import org.eclipse.elk.graph.ElkPort
import org.eclipse.elk.graph.impl.ElkPropertyToValueMapEntryImpl
import org.eclipse.elk.graph.json.text.services.ElkGraphJsonGrammarAccess
import org.eclipse.xtext.Keyword
import org.eclipse.xtext.formatting2.AbstractFormatter2
import org.eclipse.xtext.formatting2.IFormattableDocument
import org.eclipse.xtext.formatting2.IHiddenRegionFormatter
import org.eclipse.xtext.xbase.lib.Procedures.Procedure1

class ElkGraphJsonFormatter extends AbstractFormatter2 {

    static val Procedure1<? super IHiddenRegionFormatter> no_space = [noSpace]
    static val Procedure1<? super IHiddenRegionFormatter> one_space = [oneSpace]
    static val Procedure1<? super IHiddenRegionFormatter> new_line = [newLine]

    @Inject extension ElkGraphJsonGrammarAccess

    def dispatch void format(ElkNode node, extension IFormattableDocument document) {

        val compact = node.children.empty && node.ports.empty && node.labels.empty && node.properties.empty
        node.formatWithCommonGraphElementStyle(#[elkNodeAccess.commaKeyword_2_1_0, elkNodeAccess.commaKeyword_3],
            elkNodeAccess.leftCurlyBracketKeyword_1, elkNodeAccess.rightCurlyBracketKeyword_4,
            compact, 
            document
        )

        node.formatElkNodeChildren(document)
        node.formatElkNodePorts(document)
        node.formatElkGraphElementLabels(document)
        node.formatElkNodeEdges(document)
        node.formatLayoutOptions(document)
        node.formatJsonMember(document)

        node.children.forEach[it.format]
        node.ports.forEach[it.format]
        node.labels.forEach[it.format]
        node.containedEdges.forEach[it.format]
        node.properties.forEach[it.format]
    }

    def dispatch void format(ElkPort port, extension IFormattableDocument document) {
        
        val compact = port.properties.empty && port.labels.size < 2
        port.formatWithCommonGraphElementStyle(#[elkPortAccess.commaKeyword_2_0, elkPortAccess.commaKeyword_3],
            elkPortAccess.leftCurlyBracketKeyword_0, elkPortAccess.rightCurlyBracketKeyword_4, 
            compact,
            document
        )

        port.formatElkGraphElementLabels(document)
        port.formatLayoutOptions(document)
        port.formatJsonMember(document)

        port.labels.forEach[it.format]
        port.properties.forEach[it.format]
    }

    def dispatch void format(ElkLabel label, extension IFormattableDocument document) {
        
        val compact = label.properties.empty && label.labels.empty
        label.formatWithCommonGraphElementStyle(#[elkLabelAccess.commaKeyword_2_0, elkLabelAccess.commaKeyword_3],
            elkLabelAccess.leftCurlyBracketKeyword_0, elkLabelAccess.rightCurlyBracketKeyword_4, 
            compact, 
            document
        )

        label.formatElkGraphElementLabels(document)
        label.formatLayoutOptions(document)
        label.formatJsonMember(document)

        label.labels.forEach[it.format]
        label.properties.forEach[it.format]
    }

    def dispatch void format(ElkEdge edge, extension IFormattableDocument document) {
        
        val compact = edge.properties.empty && edge.sources.size < 2 && edge.targets.size < 2
        edge.formatWithCommonGraphElementStyle(#[elkEdgeAccess.commaKeyword_2_0, elkEdgeAccess.commaKeyword_3],
            elkEdgeAccess.leftCurlyBracketKeyword_0, elkEdgeAccess.rightCurlyBracketKeyword_4, 
            compact,
            document
        )

        // TODO sections
        edge.formatElkEdgeSources(document)
        edge.formatElkEdgeTargets(document)
        edge.formatElkGraphElementLabels(document)
        edge.formatLayoutOptions(document)
        edge.formatJsonMember(document)

        edge.labels.forEach[it.format]
        edge.properties.forEach[it.format]
    }

    def dispatch void format(ElkPropertyToValueMapEntryImpl entry, extension IFormattableDocument document) {
        entry.regionFor.keyword(propertyAccess.colonKeyword_1).prepend(no_space).append(one_space)
    }

    private def void formatWithCommonGraphElementStyle(ElkGraphElement element, List<Keyword> commaKeywords,
        Keyword openingKeyword, Keyword closingKeyword, boolean compact, extension IFormattableDocument document) {

        element.regionFor.keywords(':').forEach[prepend(no_space).append(one_space)]
        element.regionFor.keywords(commaKeywords).forEach [
            prepend(no_space)
            if (compact) append(one_space) else append(new_line)
        ]
        interior(
            element.regionFor.keyword(openingKeyword).prepend(no_space) => [
                if (compact) append(one_space) else append(new_line)
            ],
            element.regionFor.keyword(closingKeyword) => [
                if (compact) prepend(one_space) else prepend(new_line)
            ],
            [indent]
        )
    }

    /*
     * Formatting of container fragments.
     */
     
     private def void formatElkEdgeSources(ElkEdge edge, extension IFormattableDocument document) {
        edge.formatWithCommonContainerStyle(true,
            #[elkEdgeSourcesAccess.commaKeyword_1_1_0, elkEdgeSourcesAccess.commaKeyword_2],
            elkEdgeSourcesAccess.leftSquareBracketKeyword_0, elkEdgeSourcesAccess.rightSquareBracketKeyword_3, document)
    }

    private def void formatElkEdgeTargets(ElkEdge edge, extension IFormattableDocument document) {
        edge.formatWithCommonContainerStyle(true,
            #[elkEdgeTargetsAccess.commaKeyword_1_1_0, elkEdgeTargetsAccess.commaKeyword_2],
            elkEdgeTargetsAccess.leftSquareBracketKeyword_0, elkEdgeTargetsAccess.rightSquareBracketKeyword_3, document)
    }
     
    private def void formatElkNodeChildren(ElkNode node, extension IFormattableDocument document) {
        node.formatWithCommonContainerStyle(node.children.nullOrEmpty,
            #[elkNodeChildrenAccess.commaKeyword_1_1_0, elkNodeChildrenAccess.commaKeyword_2],
            elkNodeChildrenAccess.leftSquareBracketKeyword_0, elkNodeChildrenAccess.rightSquareBracketKeyword_3,
            document)
    }

    private def void formatElkNodePorts(ElkNode node, extension IFormattableDocument document) {
        node.formatWithCommonContainerStyle(node.ports.nullOrEmpty,
            #[elkNodePortsAccess.commaKeyword_1_1_0, elkNodePortsAccess.commaKeyword_2],
            elkNodePortsAccess.leftSquareBracketKeyword_0, elkNodePortsAccess.rightSquareBracketKeyword_3, document)
    }

    private def void formatElkNodeEdges(ElkNode node, extension IFormattableDocument document) {
        node.formatWithCommonContainerStyle(node.containedEdges.nullOrEmpty,
            #[elkNodeEdgesAccess.commaKeyword_1_1_0, elkNodeEdgesAccess.commaKeyword_2],
            elkNodeEdgesAccess.leftSquareBracketKeyword_0, elkNodeEdgesAccess.rightSquareBracketKeyword_3, document)
    }

    private def void formatElkGraphElementLabels(ElkGraphElement element, extension IFormattableDocument document) {
        element.formatWithCommonContainerStyle(element.labels.nullOrEmpty,
            #[elkGraphElementLabelsAccess.commaKeyword_1_1_0, elkGraphElementLabelsAccess.commaKeyword_2],
            elkGraphElementLabelsAccess.leftSquareBracketKeyword_0,
            elkGraphElementLabelsAccess.rightSquareBracketKeyword_3, document)
    }

    private def void formatLayoutOptions(ElkGraphElement element, extension IFormattableDocument document) {
        val compact = element.properties.size < 2
        val commas = #[elkGraphElementPropertiesAccess.commaKeyword_1_1_0,
            elkGraphElementPropertiesAccess.commaKeyword_2]
        element.formatWithCommonContainerStyle(compact, commas,
            elkGraphElementPropertiesAccess.leftCurlyBracketKeyword_0,
            elkGraphElementPropertiesAccess.rightCurlyBracketKeyword_3, document)

        // Different from the other containers
        if (!compact) {
            element.regionFor.keyword(elkGraphElementPropertiesAccess.leftCurlyBracketKeyword_0).append(new_line)
        }
        element.regionFor.keywords(commas).forEach[
            if (compact) 
                prepend(no_space).append(one_space)
            else
                append(new_line)
        ]
    }

    private def void formatWithCommonContainerStyle(ElkGraphElement element, boolean compact,
        List<Keyword> commaKeywords, Keyword openingKeyword, Keyword closingKeyword,
        extension IFormattableDocument document) {

        element.regionFor.keywords(commaKeywords).forEach[prepend(no_space).append(new_line)]
        if (compact) {
            element.regionFor.keyword(openingKeyword).append(one_space)
            element.regionFor.keyword(closingKeyword).prepend(one_space)
        } else {
            interior(
                element.regionFor.keyword(openingKeyword).prepend(one_space).append(new_line),
                element.regionFor.keyword(closingKeyword).prepend(new_line),
                [indent]
            )
        }
    }

    /*
     * TODO not working this way. There are not semantic regions found when using the 'keywords' accessor. 
     * Maybe this is due to the fact that nothing of the json members is part of the actual semantic model?
     */
    private def void formatJsonMember(ElkGraphElement element, extension IFormattableDocument document) {
        element.regionFor.keywords(jsonMemberAccess.colonKeyword_1).forEach[prepend(no_space).append(one_space)]
        element.regionFor.keywords(jsonArrayAccess.leftSquareBracketKeyword_0).forEach [
            prepend(one_space).append(new_line)
        ]
        element.regionFor.keywords(jsonArrayAccess.commaKeyword_1_1_0, jsonArrayAccess.commaKeyword_2).forEach [
            prepend(no_space).append(one_space)
        ]

        element.regionFor.keywords(jsonObjectAccess.leftCurlyBracketKeyword_0).forEach [
            prepend(one_space).append(new_line)
        ]
        element.regionFor.keywords(jsonObjectAccess.commaKeyword_1_1_0, jsonObjectAccess.commaKeyword_2).forEach [
            prepend(no_space).append(one_space)
        ]
    }

}
