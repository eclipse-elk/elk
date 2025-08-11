/*******************************************************************************
 * Copyright (c) 2017 TypeFox GmbH (http://www.typefox.io) and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.graph.text.formatting2

import com.google.inject.Inject
import org.eclipse.elk.graph.EMapPropertyHolder
import org.eclipse.elk.graph.ElkBendPoint
import org.eclipse.elk.graph.ElkEdge
import org.eclipse.elk.graph.ElkEdgeSection
import org.eclipse.elk.graph.ElkLabel
import org.eclipse.elk.graph.ElkNode
import org.eclipse.elk.graph.ElkPort
import org.eclipse.elk.graph.ElkShape
import org.eclipse.elk.graph.impl.ElkPropertyToValueMapEntryImpl
import org.eclipse.elk.graph.text.services.ElkGraphGrammarAccess
import org.eclipse.xtext.formatting2.AbstractFormatter2
import org.eclipse.xtext.formatting2.IFormattableDocument
import org.eclipse.xtext.formatting2.IHiddenRegionFormatter
import org.eclipse.xtext.xbase.lib.Procedures.Procedure1

import static org.eclipse.elk.graph.ElkGraphPackage.Literals.*

class ElkGraphFormatter extends AbstractFormatter2 {
    
    static val Procedure1<? super IHiddenRegionFormatter> no_space = [noSpace]
    static val Procedure1<? super IHiddenRegionFormatter> one_space = [oneSpace]
    static val Procedure1<? super IHiddenRegionFormatter> new_line = [newLine]
    static val Procedure1<? super IHiddenRegionFormatter> new_lines = [setNewLines(1, 1, 2)]
    
    @Inject extension ElkGraphGrammarAccess

    def dispatch void format(ElkNode node, extension IFormattableDocument document) {
        if (node.parent === null) {
            node.regionFor.keyword(rootNodeAccess.graphKeyword_1_0).append(one_space)
            node.regionFor.assignment(rootNodeAccess.identifierAssignment_1_1).append(new_lines)
        } else {
            node.regionFor.keyword(elkNodeAccess.nodeKeyword_0).append(one_space)
            interior(
                node.regionFor.keyword(elkNodeAccess.leftCurlyBracketKeyword_2_0).prepend(one_space).append(new_line),
                node.regionFor.keyword(elkNodeAccess.rightCurlyBracketKeyword_2_4),
                [indent]
            )
        }
        node.regionFor.keywords(':', ',').forEach[prepend(no_space).append(one_space)]
        node.formatShapeLayout(document)
        for (property : node.properties) {
            property.format
        }
        for (child : node.children) {
            child.append(new_lines).format
        }
        for (edge : node.containedEdges) {
            edge.append(new_lines).format
        }
        for (port : node.ports) {
            port.append(new_lines).format
        }
        for (label : node.labels) {
            label.append(new_lines).format
        }
    }
    
    def dispatch void format(ElkEdge edge, extension IFormattableDocument document) {
        edge.regionFor.keyword(elkEdgeAccess.edgeKeyword_0).append(one_space)
        edge.regionFor.keyword(elkEdgeAccess.hyphenMinusGreaterThanSignKeyword_4).prepend(one_space).append(one_space)
        edge.regionFor.keywords(':', ',').forEach[prepend(no_space).append(one_space)]
        interior(
            edge.regionFor.keyword(elkEdgeAccess.leftCurlyBracketKeyword_7_0).prepend(one_space).append(new_line),
            edge.regionFor.keyword(elkEdgeAccess.rightCurlyBracketKeyword_7_4),
            [indent]
        )
        edge.formatEdgeLayout(document)
        for (section : edge.sections) {
            section.format
            if (section.identifier !== null)
                section.append(new_lines)
        }
        for (property : edge.properties) {
            property.format
        }
        for (label : edge.labels) {
            label.append(new_lines).format
        }
    }
    
    def dispatch void format(ElkPort port, extension IFormattableDocument document) {
        port.regionFor.keyword(elkPortAccess.portKeyword_0).append(one_space)
        port.regionFor.keywords(':', ',').forEach[prepend(no_space).append(one_space)]
        interior(
            port.regionFor.keyword(elkPortAccess.leftCurlyBracketKeyword_2_0).prepend(one_space).append(new_line),
            port.regionFor.keyword(elkPortAccess.rightCurlyBracketKeyword_2_4),
            [indent]
        )
        port.formatShapeLayout(document)
        for (property : port.properties) {
            property.format
        }
        for (label : port.labels) {
            label.append(new_lines).format
        }
    }

    def dispatch void format(ElkLabel label, extension IFormattableDocument document) {
        label.regionFor.keyword(elkLabelAccess.labelKeyword_0).append(one_space)
        label.regionFor.keywords(':', ',').forEach[prepend(no_space).append(one_space)]
        interior(
            label.regionFor.keyword(elkLabelAccess.leftCurlyBracketKeyword_3_0).prepend(one_space).append(new_line),
            label.regionFor.keyword(elkLabelAccess.rightCurlyBracketKeyword_3_4),
            [indent]
        )
        label.formatShapeLayout(document)
        for (property : label.properties) {
            property.format
        }
        for (nestedLabel : label.labels) {
            nestedLabel.append(new_lines).format
        }
    }
    
    def dispatch void format(ElkEdgeSection section, extension IFormattableDocument document) {
        section.regionFor.keywords(':', ',').forEach[prepend(no_space).append(one_space)]
        if (section.identifier === null) {
            // This is a single edge section without 'section' keyword
            section.regionFor.assignment(elkSingleEdgeSectionAccess.incomingShapeAssignment_1_0_0_2).append(new_line)
            section.regionFor.assignment(elkSingleEdgeSectionAccess.outgoingShapeAssignment_1_0_1_2).append(new_line)
            section.regionFor.assignment(elkSingleEdgeSectionAccess.startYAssignment_1_0_2_4).append(new_line)
            section.regionFor.assignment(elkSingleEdgeSectionAccess.endYAssignment_1_0_3_4).append(new_line)
            section.regionFor.keywords(elkSingleEdgeSectionAccess.verticalLineKeyword_1_1_3_0).forEach[prepend(one_space).append(one_space)]
        } else {
            // The section has square brackets surrounding a list of properties
            section.regionFor.keyword(elkEdgeSectionAccess.sectionKeyword_0).append(one_space)
            section.regionFor.keyword(elkEdgeSectionAccess.hyphenMinusGreaterThanSignKeyword_2_0).prepend(one_space).append(one_space)
            if (section.hasAtLeastOneProperty) {
                interior(
                    section.regionFor.keyword(elkEdgeSectionAccess.leftSquareBracketKeyword_3).prepend(one_space).append(new_line),
                    section.regionFor.keyword(elkEdgeSectionAccess.rightSquareBracketKeyword_5),
                    [indent]
                )
                section.regionFor.assignment(elkEdgeSectionAccess.incomingShapeAssignment_4_0_0_2).append(new_line)
                section.regionFor.assignment(elkEdgeSectionAccess.outgoingShapeAssignment_4_0_1_2).append(new_line)
                section.regionFor.assignment(elkEdgeSectionAccess.startYAssignment_4_0_2_4).append(new_line)
                section.regionFor.assignment(elkEdgeSectionAccess.endYAssignment_4_0_3_4).append(new_line)
                section.regionFor.keywords(elkEdgeSectionAccess.verticalLineKeyword_4_1_3_0).forEach[prepend(one_space).append(one_space)]
            } else {
                // The section is empty
                section.regionFor.keyword(elkEdgeSectionAccess.leftSquareBracketKeyword_3).prepend(one_space)
                section.regionFor.keyword(elkEdgeSectionAccess.rightSquareBracketKeyword_5).prepend(no_space)
            }
        }
        section.bendPoints.lastOrNull.append(new_line)
        for (point : section.bendPoints) {
            point.format
        }
        for (property : section.properties) {
            property.format
        }
    }
    
    private def hasAtLeastOneProperty(ElkEdgeSection section) {
        section.eIsSet(ELK_EDGE_SECTION__INCOMING_SHAPE)
            || section.eIsSet(ELK_EDGE_SECTION__OUTGOING_SHAPE)
            || section.eIsSet(ELK_EDGE_SECTION__START_X) || section.eIsSet(ELK_EDGE_SECTION__START_Y)
            || section.eIsSet(ELK_EDGE_SECTION__END_X) || section.eIsSet(ELK_EDGE_SECTION__END_Y)
            || !section.bendPoints.empty
            || !section.properties.empty
    }
    
    def dispatch void format(ElkBendPoint bendPoint, extension IFormattableDocument document) {
        bendPoint.regionFor.keyword(elkBendPointAccess.commaKeyword_1).prepend(no_space).append(one_space)
    }
    
    def dispatch void format(ElkPropertyToValueMapEntryImpl entry, extension IFormattableDocument document) {
        entry.regionFor.keyword(propertyAccess.colonKeyword_1).prepend(no_space).append(one_space)
        val container = entry.eContainer as EMapPropertyHolder
        if (container.properties.indexOf(entry) == container.properties.size - 1)
            entry.append(new_lines)
        else
            entry.append(new_line)
    }
    
    private def void formatShapeLayout(ElkShape shape, extension IFormattableDocument document) {
        var propCount = 0
        if (shape.eIsSet(ELK_SHAPE__X) || shape.eIsSet(ELK_SHAPE__Y))
            propCount++
        if (shape.eIsSet(ELK_SHAPE__WIDTH) || shape.eIsSet(ELK_SHAPE__HEIGHT))
            propCount++
        if (propCount == 0) {
            shape.regionFor.keyword(shapeLayoutAccess.leftSquareBracketKeyword_1).prepend(one_space)
            shape.regionFor.keyword(shapeLayoutAccess.rightSquareBracketKeyword_3).prepend(no_space).append(new_lines)
        } else if (propCount == 1) {
            // Format in one line
            shape.regionFor.keyword(shapeLayoutAccess.leftSquareBracketKeyword_1).prepend(one_space)
            shape.regionFor.keyword(shapeLayoutAccess.rightSquareBracketKeyword_3).prepend(one_space).append(new_lines)
            shape.regionFor.keyword(shapeLayoutAccess.positionKeyword_2_0_0).prepend(one_space)
            shape.regionFor.keyword(shapeLayoutAccess.sizeKeyword_2_1_0).prepend(one_space)
        } else {
            // Format in multiple lines
            interior(
                shape.regionFor.keyword(shapeLayoutAccess.leftSquareBracketKeyword_1).prepend(one_space).append(new_line),
                shape.regionFor.keyword(shapeLayoutAccess.rightSquareBracketKeyword_3).append(new_lines),
                [indent]
            )
            shape.regionFor.assignment(shapeLayoutAccess.YAssignment_2_0_4).append(new_line)
            shape.regionFor.assignment(shapeLayoutAccess.heightAssignment_2_1_4).append(new_line)
        }
    }
    
    private def void formatEdgeLayout(ElkEdge edge, extension IFormattableDocument document) {
        val sectionCount = edge.sections.size
        if (sectionCount == 0) {
            edge.regionFor.keyword(edgeLayoutAccess.leftSquareBracketKeyword_1).prepend(one_space)
            edge.regionFor.keyword(edgeLayoutAccess.rightSquareBracketKeyword_3).prepend(no_space).append(new_lines)
        } else {
            interior(
                edge.regionFor.keyword(edgeLayoutAccess.leftSquareBracketKeyword_1).prepend(one_space).append(new_line),
                edge.regionFor.keyword(edgeLayoutAccess.rightSquareBracketKeyword_3).append(new_lines),
                [indent]
            )
        }
    }
    
}
