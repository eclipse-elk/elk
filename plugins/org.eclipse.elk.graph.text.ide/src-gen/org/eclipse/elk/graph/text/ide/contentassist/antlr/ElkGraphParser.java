/*******************************************************************************
 * Copyright (c) 2016 TypeFox GmbH (http://www.typefox.io) and others.
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
package org.eclipse.elk.graph.text.ide.contentassist.antlr;

import com.google.inject.Inject;
import java.util.HashMap;
import java.util.Map;
import org.eclipse.elk.graph.text.ide.contentassist.antlr.internal.InternalElkGraphParser;
import org.eclipse.elk.graph.text.services.ElkGraphGrammarAccess;
import org.eclipse.xtext.AbstractElement;
import org.eclipse.xtext.ide.editor.contentassist.antlr.AbstractContentAssistParser;

public class ElkGraphParser extends AbstractContentAssistParser {

	@Inject
	private ElkGraphGrammarAccess grammarAccess;

	private Map<AbstractElement, String> nameMappings;

	@Override
	protected InternalElkGraphParser createParser() {
		InternalElkGraphParser result = new InternalElkGraphParser(null);
		result.setGrammarAccess(grammarAccess);
		return result;
	}

	@Override
	protected String getRuleName(AbstractElement element) {
		if (nameMappings == null) {
			nameMappings = new HashMap<AbstractElement, String>() {
				private static final long serialVersionUID = 1L;
				{
					put(grammarAccess.getRootNodeAccess().getAlternatives_4(), "rule__RootNode__Alternatives_4");
					put(grammarAccess.getElkNodeAccess().getAlternatives_2_3(), "rule__ElkNode__Alternatives_2_3");
					put(grammarAccess.getEdgeLayoutAccess().getAlternatives_2(), "rule__EdgeLayout__Alternatives_2");
					put(grammarAccess.getNumberAccess().getAlternatives(), "rule__Number__Alternatives");
					put(grammarAccess.getPropertyAccess().getAlternatives_2(), "rule__Property__Alternatives_2");
					put(grammarAccess.getNumberValueAccess().getAlternatives(), "rule__NumberValue__Alternatives");
					put(grammarAccess.getBooleanValueAccess().getAlternatives(), "rule__BooleanValue__Alternatives");
					put(grammarAccess.getRootNodeAccess().getGroup(), "rule__RootNode__Group__0");
					put(grammarAccess.getRootNodeAccess().getGroup_1(), "rule__RootNode__Group_1__0");
					put(grammarAccess.getElkNodeAccess().getGroup(), "rule__ElkNode__Group__0");
					put(grammarAccess.getElkNodeAccess().getGroup_2(), "rule__ElkNode__Group_2__0");
					put(grammarAccess.getElkLabelAccess().getGroup(), "rule__ElkLabel__Group__0");
					put(grammarAccess.getElkLabelAccess().getGroup_1(), "rule__ElkLabel__Group_1__0");
					put(grammarAccess.getElkLabelAccess().getGroup_3(), "rule__ElkLabel__Group_3__0");
					put(grammarAccess.getElkPortAccess().getGroup(), "rule__ElkPort__Group__0");
					put(grammarAccess.getElkPortAccess().getGroup_2(), "rule__ElkPort__Group_2__0");
					put(grammarAccess.getShapeLayoutAccess().getGroup(), "rule__ShapeLayout__Group__0");
					put(grammarAccess.getShapeLayoutAccess().getGroup_2_0(), "rule__ShapeLayout__Group_2_0__0");
					put(grammarAccess.getShapeLayoutAccess().getGroup_2_1(), "rule__ShapeLayout__Group_2_1__0");
					put(grammarAccess.getElkEdgeAccess().getGroup(), "rule__ElkEdge__Group__0");
					put(grammarAccess.getElkEdgeAccess().getGroup_1(), "rule__ElkEdge__Group_1__0");
					put(grammarAccess.getElkEdgeAccess().getGroup_3(), "rule__ElkEdge__Group_3__0");
					put(grammarAccess.getElkEdgeAccess().getGroup_6(), "rule__ElkEdge__Group_6__0");
					put(grammarAccess.getElkEdgeAccess().getGroup_7(), "rule__ElkEdge__Group_7__0");
					put(grammarAccess.getEdgeLayoutAccess().getGroup(), "rule__EdgeLayout__Group__0");
					put(grammarAccess.getElkSingleEdgeSectionAccess().getGroup(), "rule__ElkSingleEdgeSection__Group__0");
					put(grammarAccess.getElkSingleEdgeSectionAccess().getGroup_1(), "rule__ElkSingleEdgeSection__Group_1__0");
					put(grammarAccess.getElkSingleEdgeSectionAccess().getGroup_1_0_0(), "rule__ElkSingleEdgeSection__Group_1_0_0__0");
					put(grammarAccess.getElkSingleEdgeSectionAccess().getGroup_1_0_1(), "rule__ElkSingleEdgeSection__Group_1_0_1__0");
					put(grammarAccess.getElkSingleEdgeSectionAccess().getGroup_1_0_2(), "rule__ElkSingleEdgeSection__Group_1_0_2__0");
					put(grammarAccess.getElkSingleEdgeSectionAccess().getGroup_1_0_3(), "rule__ElkSingleEdgeSection__Group_1_0_3__0");
					put(grammarAccess.getElkSingleEdgeSectionAccess().getGroup_1_1(), "rule__ElkSingleEdgeSection__Group_1_1__0");
					put(grammarAccess.getElkSingleEdgeSectionAccess().getGroup_1_1_3(), "rule__ElkSingleEdgeSection__Group_1_1_3__0");
					put(grammarAccess.getElkEdgeSectionAccess().getGroup(), "rule__ElkEdgeSection__Group__0");
					put(grammarAccess.getElkEdgeSectionAccess().getGroup_2(), "rule__ElkEdgeSection__Group_2__0");
					put(grammarAccess.getElkEdgeSectionAccess().getGroup_2_2(), "rule__ElkEdgeSection__Group_2_2__0");
					put(grammarAccess.getElkEdgeSectionAccess().getGroup_4(), "rule__ElkEdgeSection__Group_4__0");
					put(grammarAccess.getElkEdgeSectionAccess().getGroup_4_0_0(), "rule__ElkEdgeSection__Group_4_0_0__0");
					put(grammarAccess.getElkEdgeSectionAccess().getGroup_4_0_1(), "rule__ElkEdgeSection__Group_4_0_1__0");
					put(grammarAccess.getElkEdgeSectionAccess().getGroup_4_0_2(), "rule__ElkEdgeSection__Group_4_0_2__0");
					put(grammarAccess.getElkEdgeSectionAccess().getGroup_4_0_3(), "rule__ElkEdgeSection__Group_4_0_3__0");
					put(grammarAccess.getElkEdgeSectionAccess().getGroup_4_1(), "rule__ElkEdgeSection__Group_4_1__0");
					put(grammarAccess.getElkEdgeSectionAccess().getGroup_4_1_3(), "rule__ElkEdgeSection__Group_4_1_3__0");
					put(grammarAccess.getElkBendPointAccess().getGroup(), "rule__ElkBendPoint__Group__0");
					put(grammarAccess.getQualifiedIdAccess().getGroup(), "rule__QualifiedId__Group__0");
					put(grammarAccess.getQualifiedIdAccess().getGroup_1(), "rule__QualifiedId__Group_1__0");
					put(grammarAccess.getPropertyAccess().getGroup(), "rule__Property__Group__0");
					put(grammarAccess.getPropertyKeyAccess().getGroup(), "rule__PropertyKey__Group__0");
					put(grammarAccess.getPropertyKeyAccess().getGroup_1(), "rule__PropertyKey__Group_1__0");
					put(grammarAccess.getRootNodeAccess().getIdentifierAssignment_1_1(), "rule__RootNode__IdentifierAssignment_1_1");
					put(grammarAccess.getRootNodeAccess().getPropertiesAssignment_3(), "rule__RootNode__PropertiesAssignment_3");
					put(grammarAccess.getRootNodeAccess().getLabelsAssignment_4_0(), "rule__RootNode__LabelsAssignment_4_0");
					put(grammarAccess.getRootNodeAccess().getPortsAssignment_4_1(), "rule__RootNode__PortsAssignment_4_1");
					put(grammarAccess.getRootNodeAccess().getChildrenAssignment_4_2(), "rule__RootNode__ChildrenAssignment_4_2");
					put(grammarAccess.getRootNodeAccess().getContainedEdgesAssignment_4_3(), "rule__RootNode__ContainedEdgesAssignment_4_3");
					put(grammarAccess.getElkNodeAccess().getIdentifierAssignment_1(), "rule__ElkNode__IdentifierAssignment_1");
					put(grammarAccess.getElkNodeAccess().getPropertiesAssignment_2_2(), "rule__ElkNode__PropertiesAssignment_2_2");
					put(grammarAccess.getElkNodeAccess().getLabelsAssignment_2_3_0(), "rule__ElkNode__LabelsAssignment_2_3_0");
					put(grammarAccess.getElkNodeAccess().getPortsAssignment_2_3_1(), "rule__ElkNode__PortsAssignment_2_3_1");
					put(grammarAccess.getElkNodeAccess().getChildrenAssignment_2_3_2(), "rule__ElkNode__ChildrenAssignment_2_3_2");
					put(grammarAccess.getElkNodeAccess().getContainedEdgesAssignment_2_3_3(), "rule__ElkNode__ContainedEdgesAssignment_2_3_3");
					put(grammarAccess.getElkLabelAccess().getIdentifierAssignment_1_0(), "rule__ElkLabel__IdentifierAssignment_1_0");
					put(grammarAccess.getElkLabelAccess().getTextAssignment_2(), "rule__ElkLabel__TextAssignment_2");
					put(grammarAccess.getElkLabelAccess().getPropertiesAssignment_3_2(), "rule__ElkLabel__PropertiesAssignment_3_2");
					put(grammarAccess.getElkLabelAccess().getLabelsAssignment_3_3(), "rule__ElkLabel__LabelsAssignment_3_3");
					put(grammarAccess.getElkPortAccess().getIdentifierAssignment_1(), "rule__ElkPort__IdentifierAssignment_1");
					put(grammarAccess.getElkPortAccess().getPropertiesAssignment_2_2(), "rule__ElkPort__PropertiesAssignment_2_2");
					put(grammarAccess.getElkPortAccess().getLabelsAssignment_2_3(), "rule__ElkPort__LabelsAssignment_2_3");
					put(grammarAccess.getShapeLayoutAccess().getXAssignment_2_0_2(), "rule__ShapeLayout__XAssignment_2_0_2");
					put(grammarAccess.getShapeLayoutAccess().getYAssignment_2_0_4(), "rule__ShapeLayout__YAssignment_2_0_4");
					put(grammarAccess.getShapeLayoutAccess().getWidthAssignment_2_1_2(), "rule__ShapeLayout__WidthAssignment_2_1_2");
					put(grammarAccess.getShapeLayoutAccess().getHeightAssignment_2_1_4(), "rule__ShapeLayout__HeightAssignment_2_1_4");
					put(grammarAccess.getElkEdgeAccess().getIdentifierAssignment_1_0(), "rule__ElkEdge__IdentifierAssignment_1_0");
					put(grammarAccess.getElkEdgeAccess().getSourcesAssignment_2(), "rule__ElkEdge__SourcesAssignment_2");
					put(grammarAccess.getElkEdgeAccess().getSourcesAssignment_3_1(), "rule__ElkEdge__SourcesAssignment_3_1");
					put(grammarAccess.getElkEdgeAccess().getTargetsAssignment_5(), "rule__ElkEdge__TargetsAssignment_5");
					put(grammarAccess.getElkEdgeAccess().getTargetsAssignment_6_1(), "rule__ElkEdge__TargetsAssignment_6_1");
					put(grammarAccess.getElkEdgeAccess().getPropertiesAssignment_7_2(), "rule__ElkEdge__PropertiesAssignment_7_2");
					put(grammarAccess.getElkEdgeAccess().getLabelsAssignment_7_3(), "rule__ElkEdge__LabelsAssignment_7_3");
					put(grammarAccess.getEdgeLayoutAccess().getSectionsAssignment_2_0(), "rule__EdgeLayout__SectionsAssignment_2_0");
					put(grammarAccess.getEdgeLayoutAccess().getSectionsAssignment_2_1(), "rule__EdgeLayout__SectionsAssignment_2_1");
					put(grammarAccess.getElkSingleEdgeSectionAccess().getIncomingShapeAssignment_1_0_0_2(), "rule__ElkSingleEdgeSection__IncomingShapeAssignment_1_0_0_2");
					put(grammarAccess.getElkSingleEdgeSectionAccess().getOutgoingShapeAssignment_1_0_1_2(), "rule__ElkSingleEdgeSection__OutgoingShapeAssignment_1_0_1_2");
					put(grammarAccess.getElkSingleEdgeSectionAccess().getStartXAssignment_1_0_2_2(), "rule__ElkSingleEdgeSection__StartXAssignment_1_0_2_2");
					put(grammarAccess.getElkSingleEdgeSectionAccess().getStartYAssignment_1_0_2_4(), "rule__ElkSingleEdgeSection__StartYAssignment_1_0_2_4");
					put(grammarAccess.getElkSingleEdgeSectionAccess().getEndXAssignment_1_0_3_2(), "rule__ElkSingleEdgeSection__EndXAssignment_1_0_3_2");
					put(grammarAccess.getElkSingleEdgeSectionAccess().getEndYAssignment_1_0_3_4(), "rule__ElkSingleEdgeSection__EndYAssignment_1_0_3_4");
					put(grammarAccess.getElkSingleEdgeSectionAccess().getBendPointsAssignment_1_1_2(), "rule__ElkSingleEdgeSection__BendPointsAssignment_1_1_2");
					put(grammarAccess.getElkSingleEdgeSectionAccess().getBendPointsAssignment_1_1_3_1(), "rule__ElkSingleEdgeSection__BendPointsAssignment_1_1_3_1");
					put(grammarAccess.getElkSingleEdgeSectionAccess().getPropertiesAssignment_1_2(), "rule__ElkSingleEdgeSection__PropertiesAssignment_1_2");
					put(grammarAccess.getElkEdgeSectionAccess().getIdentifierAssignment_1(), "rule__ElkEdgeSection__IdentifierAssignment_1");
					put(grammarAccess.getElkEdgeSectionAccess().getOutgoingSectionsAssignment_2_1(), "rule__ElkEdgeSection__OutgoingSectionsAssignment_2_1");
					put(grammarAccess.getElkEdgeSectionAccess().getOutgoingSectionsAssignment_2_2_1(), "rule__ElkEdgeSection__OutgoingSectionsAssignment_2_2_1");
					put(grammarAccess.getElkEdgeSectionAccess().getIncomingShapeAssignment_4_0_0_2(), "rule__ElkEdgeSection__IncomingShapeAssignment_4_0_0_2");
					put(grammarAccess.getElkEdgeSectionAccess().getOutgoingShapeAssignment_4_0_1_2(), "rule__ElkEdgeSection__OutgoingShapeAssignment_4_0_1_2");
					put(grammarAccess.getElkEdgeSectionAccess().getStartXAssignment_4_0_2_2(), "rule__ElkEdgeSection__StartXAssignment_4_0_2_2");
					put(grammarAccess.getElkEdgeSectionAccess().getStartYAssignment_4_0_2_4(), "rule__ElkEdgeSection__StartYAssignment_4_0_2_4");
					put(grammarAccess.getElkEdgeSectionAccess().getEndXAssignment_4_0_3_2(), "rule__ElkEdgeSection__EndXAssignment_4_0_3_2");
					put(grammarAccess.getElkEdgeSectionAccess().getEndYAssignment_4_0_3_4(), "rule__ElkEdgeSection__EndYAssignment_4_0_3_4");
					put(grammarAccess.getElkEdgeSectionAccess().getBendPointsAssignment_4_1_2(), "rule__ElkEdgeSection__BendPointsAssignment_4_1_2");
					put(grammarAccess.getElkEdgeSectionAccess().getBendPointsAssignment_4_1_3_1(), "rule__ElkEdgeSection__BendPointsAssignment_4_1_3_1");
					put(grammarAccess.getElkEdgeSectionAccess().getPropertiesAssignment_4_2(), "rule__ElkEdgeSection__PropertiesAssignment_4_2");
					put(grammarAccess.getElkBendPointAccess().getXAssignment_0(), "rule__ElkBendPoint__XAssignment_0");
					put(grammarAccess.getElkBendPointAccess().getYAssignment_2(), "rule__ElkBendPoint__YAssignment_2");
					put(grammarAccess.getPropertyAccess().getKeyAssignment_0(), "rule__Property__KeyAssignment_0");
					put(grammarAccess.getPropertyAccess().getValueAssignment_2_0(), "rule__Property__ValueAssignment_2_0");
					put(grammarAccess.getPropertyAccess().getValueAssignment_2_1(), "rule__Property__ValueAssignment_2_1");
					put(grammarAccess.getPropertyAccess().getValueAssignment_2_2(), "rule__Property__ValueAssignment_2_2");
					put(grammarAccess.getPropertyAccess().getValueAssignment_2_3(), "rule__Property__ValueAssignment_2_3");
					put(grammarAccess.getShapeLayoutAccess().getUnorderedGroup_2(), "rule__ShapeLayout__UnorderedGroup_2");
					put(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), "rule__ElkSingleEdgeSection__UnorderedGroup_1_0");
					put(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), "rule__ElkEdgeSection__UnorderedGroup_4_0");
				}
			};
		}
		return nameMappings.get(element);
	}
			
	@Override
	protected String[] getInitialHiddenTokens() {
		return new String[] { "RULE_WS", "RULE_ML_COMMENT", "RULE_SL_COMMENT" };
	}

	public ElkGraphGrammarAccess getGrammarAccess() {
		return this.grammarAccess;
	}

	public void setGrammarAccess(ElkGraphGrammarAccess grammarAccess) {
		this.grammarAccess = grammarAccess;
	}
}
