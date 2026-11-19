/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
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
package org.eclipse.elk.graph.json.text.serializer;

import com.google.inject.Inject;
import java.util.List;
import org.eclipse.elk.graph.json.text.services.ElkGraphJsonGrammarAccess;
import org.eclipse.emf.ecore.EObject;
import org.eclipse.xtext.IGrammarAccess;
import org.eclipse.xtext.RuleCall;
import org.eclipse.xtext.nodemodel.INode;
import org.eclipse.xtext.serializer.analysis.GrammarAlias.AbstractElementAlias;
import org.eclipse.xtext.serializer.analysis.GrammarAlias.GroupAlias;
import org.eclipse.xtext.serializer.analysis.GrammarAlias.TokenAlias;
import org.eclipse.xtext.serializer.analysis.ISyntacticSequencerPDAProvider.ISynNavigable;
import org.eclipse.xtext.serializer.analysis.ISyntacticSequencerPDAProvider.ISynTransition;
import org.eclipse.xtext.serializer.sequencer.AbstractSyntacticSequencer;

@SuppressWarnings("all")
public abstract class AbstractElkGraphJsonSyntacticSequencer extends AbstractSyntacticSequencer {

	protected ElkGraphJsonGrammarAccess grammarAccess;
	protected AbstractElementAlias match_EdgeElement_ElkEdge_ElkEdgeSources___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_0_KeySourcesParserRuleCall_1_0_ColonKeyword_1_1_LeftSquareBracketKeyword_0__p;
	protected AbstractElementAlias match_EdgeElement_ElkEdge_ElkEdgeTargets___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_0_KeyTargetsParserRuleCall_2_0_ColonKeyword_2_1_LeftSquareBracketKeyword_0__p;
	protected AbstractElementAlias match_EdgeElement_ElkEdge_ElkGraphElementLabels___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_0_KeyLabelsParserRuleCall_3_0_ColonKeyword_3_1_LeftSquareBracketKeyword_0__p;
	protected AbstractElementAlias match_EdgeElement_ElkEdge_ElkGraphElementProperties___CommaKeyword_2_q_RightCurlyBracketKeyword_3_CommaKeyword_2_0_KeyLayoutOptionsParserRuleCall_4_0_ColonKeyword_4_1_LeftCurlyBracketKeyword_0__p;
	protected AbstractElementAlias match_ElkEdgeSources_CommaKeyword_2_q;
	protected AbstractElementAlias match_ElkEdgeTargets_CommaKeyword_2_q;
	protected AbstractElementAlias match_ElkEdge_CommaKeyword_3_q;
	protected AbstractElementAlias match_ElkGraphElementLabels_CommaKeyword_2_q;
	protected AbstractElementAlias match_ElkGraphElementLabels_ElkLabel_LabelElement___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_0_KeyLabelsParserRuleCall_3_0_ColonKeyword_3_1_LeftSquareBracketKeyword_0__p;
	protected AbstractElementAlias match_ElkGraphElementLabels_ElkNode_NodeElement___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_1_0_KeyLabelsParserRuleCall_4_0_ColonKeyword_4_1_LeftSquareBracketKeyword_0__p;
	protected AbstractElementAlias match_ElkGraphElementLabels_ElkPort_PortElement___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_0_KeyLabelsParserRuleCall_2_0_ColonKeyword_2_1_LeftSquareBracketKeyword_0__p;
	protected AbstractElementAlias match_ElkGraphElementProperties_CommaKeyword_2_q;
	protected AbstractElementAlias match_ElkGraphElementProperties_ElkLabel_LabelElement___CommaKeyword_2_q_RightCurlyBracketKeyword_3_CommaKeyword_2_0_KeyLayoutOptionsParserRuleCall_4_0_ColonKeyword_4_1_LeftCurlyBracketKeyword_0__p;
	protected AbstractElementAlias match_ElkGraphElementProperties_ElkNode_NodeElement___CommaKeyword_2_q_RightCurlyBracketKeyword_3_CommaKeyword_2_1_0_KeyLayoutOptionsParserRuleCall_6_0_ColonKeyword_6_1_LeftCurlyBracketKeyword_0__p;
	protected AbstractElementAlias match_ElkGraphElementProperties_ElkPort_PortElement___CommaKeyword_2_q_RightCurlyBracketKeyword_3_CommaKeyword_2_0_KeyLayoutOptionsParserRuleCall_3_0_ColonKeyword_3_1_LeftCurlyBracketKeyword_0__p;
	protected AbstractElementAlias match_ElkLabel_CommaKeyword_3_q;
	protected AbstractElementAlias match_ElkLabel_LabelElement___JsonMemberParserRuleCall_5_CommaKeyword_2_0__a;
	protected AbstractElementAlias match_ElkNodeChildren_CommaKeyword_2_q;
	protected AbstractElementAlias match_ElkNodeEdges_CommaKeyword_2_q;
	protected AbstractElementAlias match_ElkNodePorts_CommaKeyword_2_q;
	protected AbstractElementAlias match_ElkNode_CommaKeyword_3_q;
	protected AbstractElementAlias match_ElkNode_ElkNodeChildren_NodeElement___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_1_0_KeyChildrenParserRuleCall_2_0_ColonKeyword_2_1_LeftSquareBracketKeyword_0__p;
	protected AbstractElementAlias match_ElkNode_ElkNodeEdges_NodeElement___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_1_0_KeyEdgesParserRuleCall_5_0_ColonKeyword_5_1_LeftSquareBracketKeyword_0__p;
	protected AbstractElementAlias match_ElkNode_ElkNodePorts_NodeElement___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_1_0_KeyPortsParserRuleCall_3_0_ColonKeyword_3_1_LeftSquareBracketKeyword_0__p;
	protected AbstractElementAlias match_ElkPort_CommaKeyword_3_q;
	
	@Inject
	protected void init(IGrammarAccess access) {
		grammarAccess = (ElkGraphJsonGrammarAccess) access;
		match_EdgeElement_ElkEdge_ElkEdgeSources___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_0_KeySourcesParserRuleCall_1_0_ColonKeyword_1_1_LeftSquareBracketKeyword_0__p = new GroupAlias(true, false, new TokenAlias(false, true, grammarAccess.getElkEdgeSourcesAccess().getCommaKeyword_2()), new TokenAlias(false, false, grammarAccess.getElkEdgeSourcesAccess().getRightSquareBracketKeyword_3()), new TokenAlias(false, false, grammarAccess.getElkEdgeAccess().getCommaKeyword_2_0()), new TokenAlias(false, false, grammarAccess.getEdgeElementAccess().getKeySourcesParserRuleCall_1_0()), new TokenAlias(false, false, grammarAccess.getEdgeElementAccess().getColonKeyword_1_1()), new TokenAlias(false, false, grammarAccess.getElkEdgeSourcesAccess().getLeftSquareBracketKeyword_0()));
		match_EdgeElement_ElkEdge_ElkEdgeTargets___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_0_KeyTargetsParserRuleCall_2_0_ColonKeyword_2_1_LeftSquareBracketKeyword_0__p = new GroupAlias(true, false, new TokenAlias(false, true, grammarAccess.getElkEdgeTargetsAccess().getCommaKeyword_2()), new TokenAlias(false, false, grammarAccess.getElkEdgeTargetsAccess().getRightSquareBracketKeyword_3()), new TokenAlias(false, false, grammarAccess.getElkEdgeAccess().getCommaKeyword_2_0()), new TokenAlias(false, false, grammarAccess.getEdgeElementAccess().getKeyTargetsParserRuleCall_2_0()), new TokenAlias(false, false, grammarAccess.getEdgeElementAccess().getColonKeyword_2_1()), new TokenAlias(false, false, grammarAccess.getElkEdgeTargetsAccess().getLeftSquareBracketKeyword_0()));
		match_EdgeElement_ElkEdge_ElkGraphElementLabels___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_0_KeyLabelsParserRuleCall_3_0_ColonKeyword_3_1_LeftSquareBracketKeyword_0__p = new GroupAlias(true, false, new TokenAlias(false, true, grammarAccess.getElkGraphElementLabelsAccess().getCommaKeyword_2()), new TokenAlias(false, false, grammarAccess.getElkGraphElementLabelsAccess().getRightSquareBracketKeyword_3()), new TokenAlias(false, false, grammarAccess.getElkEdgeAccess().getCommaKeyword_2_0()), new TokenAlias(false, false, grammarAccess.getEdgeElementAccess().getKeyLabelsParserRuleCall_3_0()), new TokenAlias(false, false, grammarAccess.getEdgeElementAccess().getColonKeyword_3_1()), new TokenAlias(false, false, grammarAccess.getElkGraphElementLabelsAccess().getLeftSquareBracketKeyword_0()));
		match_EdgeElement_ElkEdge_ElkGraphElementProperties___CommaKeyword_2_q_RightCurlyBracketKeyword_3_CommaKeyword_2_0_KeyLayoutOptionsParserRuleCall_4_0_ColonKeyword_4_1_LeftCurlyBracketKeyword_0__p = new GroupAlias(true, false, new TokenAlias(false, true, grammarAccess.getElkGraphElementPropertiesAccess().getCommaKeyword_2()), new TokenAlias(false, false, grammarAccess.getElkGraphElementPropertiesAccess().getRightCurlyBracketKeyword_3()), new TokenAlias(false, false, grammarAccess.getElkEdgeAccess().getCommaKeyword_2_0()), new TokenAlias(false, false, grammarAccess.getEdgeElementAccess().getKeyLayoutOptionsParserRuleCall_4_0()), new TokenAlias(false, false, grammarAccess.getEdgeElementAccess().getColonKeyword_4_1()), new TokenAlias(false, false, grammarAccess.getElkGraphElementPropertiesAccess().getLeftCurlyBracketKeyword_0()));
		match_ElkEdgeSources_CommaKeyword_2_q = new TokenAlias(false, true, grammarAccess.getElkEdgeSourcesAccess().getCommaKeyword_2());
		match_ElkEdgeTargets_CommaKeyword_2_q = new TokenAlias(false, true, grammarAccess.getElkEdgeTargetsAccess().getCommaKeyword_2());
		match_ElkEdge_CommaKeyword_3_q = new TokenAlias(false, true, grammarAccess.getElkEdgeAccess().getCommaKeyword_3());
		match_ElkGraphElementLabels_CommaKeyword_2_q = new TokenAlias(false, true, grammarAccess.getElkGraphElementLabelsAccess().getCommaKeyword_2());
		match_ElkGraphElementLabels_ElkLabel_LabelElement___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_0_KeyLabelsParserRuleCall_3_0_ColonKeyword_3_1_LeftSquareBracketKeyword_0__p = new GroupAlias(true, false, new TokenAlias(false, true, grammarAccess.getElkGraphElementLabelsAccess().getCommaKeyword_2()), new TokenAlias(false, false, grammarAccess.getElkGraphElementLabelsAccess().getRightSquareBracketKeyword_3()), new TokenAlias(false, false, grammarAccess.getElkLabelAccess().getCommaKeyword_2_0()), new TokenAlias(false, false, grammarAccess.getLabelElementAccess().getKeyLabelsParserRuleCall_3_0()), new TokenAlias(false, false, grammarAccess.getLabelElementAccess().getColonKeyword_3_1()), new TokenAlias(false, false, grammarAccess.getElkGraphElementLabelsAccess().getLeftSquareBracketKeyword_0()));
		match_ElkGraphElementLabels_ElkNode_NodeElement___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_1_0_KeyLabelsParserRuleCall_4_0_ColonKeyword_4_1_LeftSquareBracketKeyword_0__p = new GroupAlias(true, false, new TokenAlias(false, true, grammarAccess.getElkGraphElementLabelsAccess().getCommaKeyword_2()), new TokenAlias(false, false, grammarAccess.getElkGraphElementLabelsAccess().getRightSquareBracketKeyword_3()), new TokenAlias(false, false, grammarAccess.getElkNodeAccess().getCommaKeyword_2_1_0()), new TokenAlias(false, false, grammarAccess.getNodeElementAccess().getKeyLabelsParserRuleCall_4_0()), new TokenAlias(false, false, grammarAccess.getNodeElementAccess().getColonKeyword_4_1()), new TokenAlias(false, false, grammarAccess.getElkGraphElementLabelsAccess().getLeftSquareBracketKeyword_0()));
		match_ElkGraphElementLabels_ElkPort_PortElement___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_0_KeyLabelsParserRuleCall_2_0_ColonKeyword_2_1_LeftSquareBracketKeyword_0__p = new GroupAlias(true, false, new TokenAlias(false, true, grammarAccess.getElkGraphElementLabelsAccess().getCommaKeyword_2()), new TokenAlias(false, false, grammarAccess.getElkGraphElementLabelsAccess().getRightSquareBracketKeyword_3()), new TokenAlias(false, false, grammarAccess.getElkPortAccess().getCommaKeyword_2_0()), new TokenAlias(false, false, grammarAccess.getPortElementAccess().getKeyLabelsParserRuleCall_2_0()), new TokenAlias(false, false, grammarAccess.getPortElementAccess().getColonKeyword_2_1()), new TokenAlias(false, false, grammarAccess.getElkGraphElementLabelsAccess().getLeftSquareBracketKeyword_0()));
		match_ElkGraphElementProperties_CommaKeyword_2_q = new TokenAlias(false, true, grammarAccess.getElkGraphElementPropertiesAccess().getCommaKeyword_2());
		match_ElkGraphElementProperties_ElkLabel_LabelElement___CommaKeyword_2_q_RightCurlyBracketKeyword_3_CommaKeyword_2_0_KeyLayoutOptionsParserRuleCall_4_0_ColonKeyword_4_1_LeftCurlyBracketKeyword_0__p = new GroupAlias(true, false, new TokenAlias(false, true, grammarAccess.getElkGraphElementPropertiesAccess().getCommaKeyword_2()), new TokenAlias(false, false, grammarAccess.getElkGraphElementPropertiesAccess().getRightCurlyBracketKeyword_3()), new TokenAlias(false, false, grammarAccess.getElkLabelAccess().getCommaKeyword_2_0()), new TokenAlias(false, false, grammarAccess.getLabelElementAccess().getKeyLayoutOptionsParserRuleCall_4_0()), new TokenAlias(false, false, grammarAccess.getLabelElementAccess().getColonKeyword_4_1()), new TokenAlias(false, false, grammarAccess.getElkGraphElementPropertiesAccess().getLeftCurlyBracketKeyword_0()));
		match_ElkGraphElementProperties_ElkNode_NodeElement___CommaKeyword_2_q_RightCurlyBracketKeyword_3_CommaKeyword_2_1_0_KeyLayoutOptionsParserRuleCall_6_0_ColonKeyword_6_1_LeftCurlyBracketKeyword_0__p = new GroupAlias(true, false, new TokenAlias(false, true, grammarAccess.getElkGraphElementPropertiesAccess().getCommaKeyword_2()), new TokenAlias(false, false, grammarAccess.getElkGraphElementPropertiesAccess().getRightCurlyBracketKeyword_3()), new TokenAlias(false, false, grammarAccess.getElkNodeAccess().getCommaKeyword_2_1_0()), new TokenAlias(false, false, grammarAccess.getNodeElementAccess().getKeyLayoutOptionsParserRuleCall_6_0()), new TokenAlias(false, false, grammarAccess.getNodeElementAccess().getColonKeyword_6_1()), new TokenAlias(false, false, grammarAccess.getElkGraphElementPropertiesAccess().getLeftCurlyBracketKeyword_0()));
		match_ElkGraphElementProperties_ElkPort_PortElement___CommaKeyword_2_q_RightCurlyBracketKeyword_3_CommaKeyword_2_0_KeyLayoutOptionsParserRuleCall_3_0_ColonKeyword_3_1_LeftCurlyBracketKeyword_0__p = new GroupAlias(true, false, new TokenAlias(false, true, grammarAccess.getElkGraphElementPropertiesAccess().getCommaKeyword_2()), new TokenAlias(false, false, grammarAccess.getElkGraphElementPropertiesAccess().getRightCurlyBracketKeyword_3()), new TokenAlias(false, false, grammarAccess.getElkPortAccess().getCommaKeyword_2_0()), new TokenAlias(false, false, grammarAccess.getPortElementAccess().getKeyLayoutOptionsParserRuleCall_3_0()), new TokenAlias(false, false, grammarAccess.getPortElementAccess().getColonKeyword_3_1()), new TokenAlias(false, false, grammarAccess.getElkGraphElementPropertiesAccess().getLeftCurlyBracketKeyword_0()));
		match_ElkLabel_CommaKeyword_3_q = new TokenAlias(false, true, grammarAccess.getElkLabelAccess().getCommaKeyword_3());
		match_ElkLabel_LabelElement___JsonMemberParserRuleCall_5_CommaKeyword_2_0__a = new GroupAlias(true, true, new TokenAlias(false, false, grammarAccess.getLabelElementAccess().getJsonMemberParserRuleCall_5()), new TokenAlias(false, false, grammarAccess.getElkLabelAccess().getCommaKeyword_2_0()));
		match_ElkNodeChildren_CommaKeyword_2_q = new TokenAlias(false, true, grammarAccess.getElkNodeChildrenAccess().getCommaKeyword_2());
		match_ElkNodeEdges_CommaKeyword_2_q = new TokenAlias(false, true, grammarAccess.getElkNodeEdgesAccess().getCommaKeyword_2());
		match_ElkNodePorts_CommaKeyword_2_q = new TokenAlias(false, true, grammarAccess.getElkNodePortsAccess().getCommaKeyword_2());
		match_ElkNode_CommaKeyword_3_q = new TokenAlias(false, true, grammarAccess.getElkNodeAccess().getCommaKeyword_3());
		match_ElkNode_ElkNodeChildren_NodeElement___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_1_0_KeyChildrenParserRuleCall_2_0_ColonKeyword_2_1_LeftSquareBracketKeyword_0__p = new GroupAlias(true, false, new TokenAlias(false, true, grammarAccess.getElkNodeChildrenAccess().getCommaKeyword_2()), new TokenAlias(false, false, grammarAccess.getElkNodeChildrenAccess().getRightSquareBracketKeyword_3()), new TokenAlias(false, false, grammarAccess.getElkNodeAccess().getCommaKeyword_2_1_0()), new TokenAlias(false, false, grammarAccess.getNodeElementAccess().getKeyChildrenParserRuleCall_2_0()), new TokenAlias(false, false, grammarAccess.getNodeElementAccess().getColonKeyword_2_1()), new TokenAlias(false, false, grammarAccess.getElkNodeChildrenAccess().getLeftSquareBracketKeyword_0()));
		match_ElkNode_ElkNodeEdges_NodeElement___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_1_0_KeyEdgesParserRuleCall_5_0_ColonKeyword_5_1_LeftSquareBracketKeyword_0__p = new GroupAlias(true, false, new TokenAlias(false, true, grammarAccess.getElkNodeEdgesAccess().getCommaKeyword_2()), new TokenAlias(false, false, grammarAccess.getElkNodeEdgesAccess().getRightSquareBracketKeyword_3()), new TokenAlias(false, false, grammarAccess.getElkNodeAccess().getCommaKeyword_2_1_0()), new TokenAlias(false, false, grammarAccess.getNodeElementAccess().getKeyEdgesParserRuleCall_5_0()), new TokenAlias(false, false, grammarAccess.getNodeElementAccess().getColonKeyword_5_1()), new TokenAlias(false, false, grammarAccess.getElkNodeEdgesAccess().getLeftSquareBracketKeyword_0()));
		match_ElkNode_ElkNodePorts_NodeElement___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_1_0_KeyPortsParserRuleCall_3_0_ColonKeyword_3_1_LeftSquareBracketKeyword_0__p = new GroupAlias(true, false, new TokenAlias(false, true, grammarAccess.getElkNodePortsAccess().getCommaKeyword_2()), new TokenAlias(false, false, grammarAccess.getElkNodePortsAccess().getRightSquareBracketKeyword_3()), new TokenAlias(false, false, grammarAccess.getElkNodeAccess().getCommaKeyword_2_1_0()), new TokenAlias(false, false, grammarAccess.getNodeElementAccess().getKeyPortsParserRuleCall_3_0()), new TokenAlias(false, false, grammarAccess.getNodeElementAccess().getColonKeyword_3_1()), new TokenAlias(false, false, grammarAccess.getElkNodePortsAccess().getLeftSquareBracketKeyword_0()));
		match_ElkPort_CommaKeyword_3_q = new TokenAlias(false, true, grammarAccess.getElkPortAccess().getCommaKeyword_3());
	}
	
	@Override
	protected String getUnassignedRuleCallToken(EObject semanticObject, RuleCall ruleCall, INode node) {
		if (ruleCall.getRule() == grammarAccess.getJsonMemberRule())
			return getJsonMemberToken(semanticObject, ruleCall, node);
		else if (ruleCall.getRule() == grammarAccess.getKeyChildrenRule())
			return getKeyChildrenToken(semanticObject, ruleCall, node);
		else if (ruleCall.getRule() == grammarAccess.getKeyEdgesRule())
			return getKeyEdgesToken(semanticObject, ruleCall, node);
		else if (ruleCall.getRule() == grammarAccess.getKeyHeightRule())
			return getKeyHeightToken(semanticObject, ruleCall, node);
		else if (ruleCall.getRule() == grammarAccess.getKeyIdRule())
			return getKeyIdToken(semanticObject, ruleCall, node);
		else if (ruleCall.getRule() == grammarAccess.getKeyLabelsRule())
			return getKeyLabelsToken(semanticObject, ruleCall, node);
		else if (ruleCall.getRule() == grammarAccess.getKeyLayoutOptionsRule())
			return getKeyLayoutOptionsToken(semanticObject, ruleCall, node);
		else if (ruleCall.getRule() == grammarAccess.getKeyPortsRule())
			return getKeyPortsToken(semanticObject, ruleCall, node);
		else if (ruleCall.getRule() == grammarAccess.getKeySourcesRule())
			return getKeySourcesToken(semanticObject, ruleCall, node);
		else if (ruleCall.getRule() == grammarAccess.getKeyTargetsRule())
			return getKeyTargetsToken(semanticObject, ruleCall, node);
		else if (ruleCall.getRule() == grammarAccess.getKeyTextRule())
			return getKeyTextToken(semanticObject, ruleCall, node);
		else if (ruleCall.getRule() == grammarAccess.getKeyWidthRule())
			return getKeyWidthToken(semanticObject, ruleCall, node);
		else if (ruleCall.getRule() == grammarAccess.getKeyXRule())
			return getKeyXToken(semanticObject, ruleCall, node);
		else if (ruleCall.getRule() == grammarAccess.getKeyYRule())
			return getKeyYToken(semanticObject, ruleCall, node);
		return "";
	}
	
	/**
	 * JsonMember: (STRING | ID) ':' JsonValue;
	 */
	protected String getJsonMemberToken(EObject semanticObject, RuleCall ruleCall, INode node) {
		if (node != null)
			return getTokenText(node);
		return "\"\" :";
	}
	
	/**
	 * KeyChildren: '"children"' | "'children'" | 'children';
	 */
	protected String getKeyChildrenToken(EObject semanticObject, RuleCall ruleCall, INode node) {
		if (node != null)
			return getTokenText(node);
		return "\"children\"";
	}
	
	/**
	 * KeyEdges: '"edges"' | "'edges'" | 'edges';
	 */
	protected String getKeyEdgesToken(EObject semanticObject, RuleCall ruleCall, INode node) {
		if (node != null)
			return getTokenText(node);
		return "\"edges\"";
	}
	
	/**
	 * KeyHeight: '"height"' | "'height'" | 'height';
	 */
	protected String getKeyHeightToken(EObject semanticObject, RuleCall ruleCall, INode node) {
		if (node != null)
			return getTokenText(node);
		return "\"height\"";
	}
	
	/**
	 * KeyId: '"id"' | "'id'" | 'id';
	 */
	protected String getKeyIdToken(EObject semanticObject, RuleCall ruleCall, INode node) {
		if (node != null)
			return getTokenText(node);
		return "\"id\"";
	}
	
	/**
	 * KeyLabels: '"labels"' | "'labels'" | 'labels';
	 */
	protected String getKeyLabelsToken(EObject semanticObject, RuleCall ruleCall, INode node) {
		if (node != null)
			return getTokenText(node);
		return "\"labels\"";
	}
	
	/**
	 * KeyLayoutOptions:
	 *     '"layoutOptions"' | "'layoutOptions'" | 'layoutOptions'
	 *     | '"properties"' | "'properties'" | 'properties'
	 * ;
	 */
	protected String getKeyLayoutOptionsToken(EObject semanticObject, RuleCall ruleCall, INode node) {
		if (node != null)
			return getTokenText(node);
		return "\"layoutOptions\"";
	}
	
	/**
	 * KeyPorts: '"ports"' | "'ports'" | 'ports';
	 */
	protected String getKeyPortsToken(EObject semanticObject, RuleCall ruleCall, INode node) {
		if (node != null)
			return getTokenText(node);
		return "\"ports\"";
	}
	
	/**
	 * KeySources: '"sources"' | "'sources'" | 'sources';
	 */
	protected String getKeySourcesToken(EObject semanticObject, RuleCall ruleCall, INode node) {
		if (node != null)
			return getTokenText(node);
		return "\"sources\"";
	}
	
	/**
	 * KeyTargets: '"targets"' | "'targets'" | 'targets';
	 */
	protected String getKeyTargetsToken(EObject semanticObject, RuleCall ruleCall, INode node) {
		if (node != null)
			return getTokenText(node);
		return "\"targets\"";
	}
	
	/**
	 * KeyText: '"text"' | "'text'" | 'text';
	 */
	protected String getKeyTextToken(EObject semanticObject, RuleCall ruleCall, INode node) {
		if (node != null)
			return getTokenText(node);
		return "\"text\"";
	}
	
	/**
	 * KeyWidth: '"width"' | "'width'" | 'width';
	 */
	protected String getKeyWidthToken(EObject semanticObject, RuleCall ruleCall, INode node) {
		if (node != null)
			return getTokenText(node);
		return "\"width\"";
	}
	
	/**
	 * KeyX: '"x"' | "'x'" | 'x';
	 */
	protected String getKeyXToken(EObject semanticObject, RuleCall ruleCall, INode node) {
		if (node != null)
			return getTokenText(node);
		return "\"x\"";
	}
	
	/**
	 * KeyY: '"y"' | "'y'" | 'y';
	 */
	protected String getKeyYToken(EObject semanticObject, RuleCall ruleCall, INode node) {
		if (node != null)
			return getTokenText(node);
		return "\"y\"";
	}
	
	@Override
	protected void emitUnassignedTokens(EObject semanticObject, ISynTransition transition, INode fromNode, INode toNode) {
		if (transition.getAmbiguousSyntaxes().isEmpty()) return;
		List<INode> transitionNodes = collectNodes(fromNode, toNode);
		for (AbstractElementAlias syntax : transition.getAmbiguousSyntaxes()) {
			List<INode> syntaxNodes = getNodesFor(transitionNodes, syntax);
			if (match_EdgeElement_ElkEdge_ElkEdgeSources___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_0_KeySourcesParserRuleCall_1_0_ColonKeyword_1_1_LeftSquareBracketKeyword_0__p.equals(syntax))
				emit_EdgeElement_ElkEdge_ElkEdgeSources___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_0_KeySourcesParserRuleCall_1_0_ColonKeyword_1_1_LeftSquareBracketKeyword_0__p(semanticObject, getLastNavigableState(), syntaxNodes);
			else if (match_EdgeElement_ElkEdge_ElkEdgeTargets___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_0_KeyTargetsParserRuleCall_2_0_ColonKeyword_2_1_LeftSquareBracketKeyword_0__p.equals(syntax))
				emit_EdgeElement_ElkEdge_ElkEdgeTargets___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_0_KeyTargetsParserRuleCall_2_0_ColonKeyword_2_1_LeftSquareBracketKeyword_0__p(semanticObject, getLastNavigableState(), syntaxNodes);
			else if (match_EdgeElement_ElkEdge_ElkGraphElementLabels___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_0_KeyLabelsParserRuleCall_3_0_ColonKeyword_3_1_LeftSquareBracketKeyword_0__p.equals(syntax))
				emit_EdgeElement_ElkEdge_ElkGraphElementLabels___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_0_KeyLabelsParserRuleCall_3_0_ColonKeyword_3_1_LeftSquareBracketKeyword_0__p(semanticObject, getLastNavigableState(), syntaxNodes);
			else if (match_EdgeElement_ElkEdge_ElkGraphElementProperties___CommaKeyword_2_q_RightCurlyBracketKeyword_3_CommaKeyword_2_0_KeyLayoutOptionsParserRuleCall_4_0_ColonKeyword_4_1_LeftCurlyBracketKeyword_0__p.equals(syntax))
				emit_EdgeElement_ElkEdge_ElkGraphElementProperties___CommaKeyword_2_q_RightCurlyBracketKeyword_3_CommaKeyword_2_0_KeyLayoutOptionsParserRuleCall_4_0_ColonKeyword_4_1_LeftCurlyBracketKeyword_0__p(semanticObject, getLastNavigableState(), syntaxNodes);
			else if (match_ElkEdgeSources_CommaKeyword_2_q.equals(syntax))
				emit_ElkEdgeSources_CommaKeyword_2_q(semanticObject, getLastNavigableState(), syntaxNodes);
			else if (match_ElkEdgeTargets_CommaKeyword_2_q.equals(syntax))
				emit_ElkEdgeTargets_CommaKeyword_2_q(semanticObject, getLastNavigableState(), syntaxNodes);
			else if (match_ElkEdge_CommaKeyword_3_q.equals(syntax))
				emit_ElkEdge_CommaKeyword_3_q(semanticObject, getLastNavigableState(), syntaxNodes);
			else if (match_ElkGraphElementLabels_CommaKeyword_2_q.equals(syntax))
				emit_ElkGraphElementLabels_CommaKeyword_2_q(semanticObject, getLastNavigableState(), syntaxNodes);
			else if (match_ElkGraphElementLabels_ElkLabel_LabelElement___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_0_KeyLabelsParserRuleCall_3_0_ColonKeyword_3_1_LeftSquareBracketKeyword_0__p.equals(syntax))
				emit_ElkGraphElementLabels_ElkLabel_LabelElement___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_0_KeyLabelsParserRuleCall_3_0_ColonKeyword_3_1_LeftSquareBracketKeyword_0__p(semanticObject, getLastNavigableState(), syntaxNodes);
			else if (match_ElkGraphElementLabels_ElkNode_NodeElement___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_1_0_KeyLabelsParserRuleCall_4_0_ColonKeyword_4_1_LeftSquareBracketKeyword_0__p.equals(syntax))
				emit_ElkGraphElementLabels_ElkNode_NodeElement___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_1_0_KeyLabelsParserRuleCall_4_0_ColonKeyword_4_1_LeftSquareBracketKeyword_0__p(semanticObject, getLastNavigableState(), syntaxNodes);
			else if (match_ElkGraphElementLabels_ElkPort_PortElement___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_0_KeyLabelsParserRuleCall_2_0_ColonKeyword_2_1_LeftSquareBracketKeyword_0__p.equals(syntax))
				emit_ElkGraphElementLabels_ElkPort_PortElement___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_0_KeyLabelsParserRuleCall_2_0_ColonKeyword_2_1_LeftSquareBracketKeyword_0__p(semanticObject, getLastNavigableState(), syntaxNodes);
			else if (match_ElkGraphElementProperties_CommaKeyword_2_q.equals(syntax))
				emit_ElkGraphElementProperties_CommaKeyword_2_q(semanticObject, getLastNavigableState(), syntaxNodes);
			else if (match_ElkGraphElementProperties_ElkLabel_LabelElement___CommaKeyword_2_q_RightCurlyBracketKeyword_3_CommaKeyword_2_0_KeyLayoutOptionsParserRuleCall_4_0_ColonKeyword_4_1_LeftCurlyBracketKeyword_0__p.equals(syntax))
				emit_ElkGraphElementProperties_ElkLabel_LabelElement___CommaKeyword_2_q_RightCurlyBracketKeyword_3_CommaKeyword_2_0_KeyLayoutOptionsParserRuleCall_4_0_ColonKeyword_4_1_LeftCurlyBracketKeyword_0__p(semanticObject, getLastNavigableState(), syntaxNodes);
			else if (match_ElkGraphElementProperties_ElkNode_NodeElement___CommaKeyword_2_q_RightCurlyBracketKeyword_3_CommaKeyword_2_1_0_KeyLayoutOptionsParserRuleCall_6_0_ColonKeyword_6_1_LeftCurlyBracketKeyword_0__p.equals(syntax))
				emit_ElkGraphElementProperties_ElkNode_NodeElement___CommaKeyword_2_q_RightCurlyBracketKeyword_3_CommaKeyword_2_1_0_KeyLayoutOptionsParserRuleCall_6_0_ColonKeyword_6_1_LeftCurlyBracketKeyword_0__p(semanticObject, getLastNavigableState(), syntaxNodes);
			else if (match_ElkGraphElementProperties_ElkPort_PortElement___CommaKeyword_2_q_RightCurlyBracketKeyword_3_CommaKeyword_2_0_KeyLayoutOptionsParserRuleCall_3_0_ColonKeyword_3_1_LeftCurlyBracketKeyword_0__p.equals(syntax))
				emit_ElkGraphElementProperties_ElkPort_PortElement___CommaKeyword_2_q_RightCurlyBracketKeyword_3_CommaKeyword_2_0_KeyLayoutOptionsParserRuleCall_3_0_ColonKeyword_3_1_LeftCurlyBracketKeyword_0__p(semanticObject, getLastNavigableState(), syntaxNodes);
			else if (match_ElkLabel_CommaKeyword_3_q.equals(syntax))
				emit_ElkLabel_CommaKeyword_3_q(semanticObject, getLastNavigableState(), syntaxNodes);
			else if (match_ElkLabel_LabelElement___JsonMemberParserRuleCall_5_CommaKeyword_2_0__a.equals(syntax))
				emit_ElkLabel_LabelElement___JsonMemberParserRuleCall_5_CommaKeyword_2_0__a(semanticObject, getLastNavigableState(), syntaxNodes);
			else if (match_ElkNodeChildren_CommaKeyword_2_q.equals(syntax))
				emit_ElkNodeChildren_CommaKeyword_2_q(semanticObject, getLastNavigableState(), syntaxNodes);
			else if (match_ElkNodeEdges_CommaKeyword_2_q.equals(syntax))
				emit_ElkNodeEdges_CommaKeyword_2_q(semanticObject, getLastNavigableState(), syntaxNodes);
			else if (match_ElkNodePorts_CommaKeyword_2_q.equals(syntax))
				emit_ElkNodePorts_CommaKeyword_2_q(semanticObject, getLastNavigableState(), syntaxNodes);
			else if (match_ElkNode_CommaKeyword_3_q.equals(syntax))
				emit_ElkNode_CommaKeyword_3_q(semanticObject, getLastNavigableState(), syntaxNodes);
			else if (match_ElkNode_ElkNodeChildren_NodeElement___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_1_0_KeyChildrenParserRuleCall_2_0_ColonKeyword_2_1_LeftSquareBracketKeyword_0__p.equals(syntax))
				emit_ElkNode_ElkNodeChildren_NodeElement___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_1_0_KeyChildrenParserRuleCall_2_0_ColonKeyword_2_1_LeftSquareBracketKeyword_0__p(semanticObject, getLastNavigableState(), syntaxNodes);
			else if (match_ElkNode_ElkNodeEdges_NodeElement___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_1_0_KeyEdgesParserRuleCall_5_0_ColonKeyword_5_1_LeftSquareBracketKeyword_0__p.equals(syntax))
				emit_ElkNode_ElkNodeEdges_NodeElement___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_1_0_KeyEdgesParserRuleCall_5_0_ColonKeyword_5_1_LeftSquareBracketKeyword_0__p(semanticObject, getLastNavigableState(), syntaxNodes);
			else if (match_ElkNode_ElkNodePorts_NodeElement___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_1_0_KeyPortsParserRuleCall_3_0_ColonKeyword_3_1_LeftSquareBracketKeyword_0__p.equals(syntax))
				emit_ElkNode_ElkNodePorts_NodeElement___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_1_0_KeyPortsParserRuleCall_3_0_ColonKeyword_3_1_LeftSquareBracketKeyword_0__p(semanticObject, getLastNavigableState(), syntaxNodes);
			else if (match_ElkPort_CommaKeyword_3_q.equals(syntax))
				emit_ElkPort_CommaKeyword_3_q(semanticObject, getLastNavigableState(), syntaxNodes);
			else acceptNodes(getLastNavigableState(), syntaxNodes);
		}
	}

	/**
	 * Ambiguous syntax:
	 *     (
	  *         ','? 
	  *         ']' 
	  *         ',' 
	  *         KeySources 
	  *         ':' 
	  *         '['
	  *     )+
	 *
	 * This ambiguous syntax occurs at:
	 *     sources+=[ElkConnectableShape|STRING] (ambiguity) sources+=[ElkConnectableShape|STRING]
	 */
	protected void emit_EdgeElement_ElkEdge_ElkEdgeSources___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_0_KeySourcesParserRuleCall_1_0_ColonKeyword_1_1_LeftSquareBracketKeyword_0__p(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
	/**
	 * Ambiguous syntax:
	 *     (
	  *         ','? 
	  *         ']' 
	  *         ',' 
	  *         KeyTargets 
	  *         ':' 
	  *         '['
	  *     )+
	 *
	 * This ambiguous syntax occurs at:
	 *     targets+=[ElkConnectableShape|STRING] (ambiguity) targets+=[ElkConnectableShape|STRING]
	 */
	protected void emit_EdgeElement_ElkEdge_ElkEdgeTargets___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_0_KeyTargetsParserRuleCall_2_0_ColonKeyword_2_1_LeftSquareBracketKeyword_0__p(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
	/**
	 * Ambiguous syntax:
	 *     (
	  *         ','? 
	  *         ']' 
	  *         ',' 
	  *         KeyLabels 
	  *         ':' 
	  *         '['
	  *     )+
	 *
	 * This ambiguous syntax occurs at:
	 *     labels+=ElkLabel (ambiguity) labels+=ElkLabel
	 */
	protected void emit_EdgeElement_ElkEdge_ElkGraphElementLabels___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_0_KeyLabelsParserRuleCall_3_0_ColonKeyword_3_1_LeftSquareBracketKeyword_0__p(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
	/**
	 * Ambiguous syntax:
	 *     (
	  *         ','? 
	  *         '}' 
	  *         ',' 
	  *         KeyLayoutOptions 
	  *         ':' 
	  *         '{'
	  *     )+
	 *
	 * This ambiguous syntax occurs at:
	 *     properties+=Property (ambiguity) properties+=Property
	 */
	protected void emit_EdgeElement_ElkEdge_ElkGraphElementProperties___CommaKeyword_2_q_RightCurlyBracketKeyword_3_CommaKeyword_2_0_KeyLayoutOptionsParserRuleCall_4_0_ColonKeyword_4_1_LeftCurlyBracketKeyword_0__p(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
	/**
	 * Ambiguous syntax:
	 *     ','?
	 *
	 * This ambiguous syntax occurs at:
	 *     sources+=[ElkConnectableShape|STRING] (ambiguity) ']' ',' KeyId ':' identifier=STRING
	 *     sources+=[ElkConnectableShape|STRING] (ambiguity) ']' ',' KeyLabels ':' '[' labels+=ElkLabel
	 *     sources+=[ElkConnectableShape|STRING] (ambiguity) ']' ',' KeyLayoutOptions ':' '{' properties+=Property
	 *     sources+=[ElkConnectableShape|STRING] (ambiguity) ']' ',' KeySources ':' '[' sources+=[ElkConnectableShape|STRING]
	 *     sources+=[ElkConnectableShape|STRING] (ambiguity) ']' ',' KeyTargets ':' '[' targets+=[ElkConnectableShape|STRING]
	 *     sources+=[ElkConnectableShape|STRING] (ambiguity) ']' ','? '}' (rule end)
	 */
	protected void emit_ElkEdgeSources_CommaKeyword_2_q(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
	/**
	 * Ambiguous syntax:
	 *     ','?
	 *
	 * This ambiguous syntax occurs at:
	 *     targets+=[ElkConnectableShape|STRING] (ambiguity) ']' ',' KeyId ':' identifier=STRING
	 *     targets+=[ElkConnectableShape|STRING] (ambiguity) ']' ',' KeyLabels ':' '[' labels+=ElkLabel
	 *     targets+=[ElkConnectableShape|STRING] (ambiguity) ']' ',' KeyLayoutOptions ':' '{' properties+=Property
	 *     targets+=[ElkConnectableShape|STRING] (ambiguity) ']' ',' KeySources ':' '[' sources+=[ElkConnectableShape|STRING]
	 *     targets+=[ElkConnectableShape|STRING] (ambiguity) ']' ',' KeyTargets ':' '[' targets+=[ElkConnectableShape|STRING]
	 *     targets+=[ElkConnectableShape|STRING] (ambiguity) ']' ','? '}' (rule end)
	 */
	protected void emit_ElkEdgeTargets_CommaKeyword_2_q(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
	/**
	 * Ambiguous syntax:
	 *     ','?
	 *
	 * This ambiguous syntax occurs at:
	 *     (rule start) '{' JsonMember (ambiguity) '}' (rule start)
	 *     identifier=STRING (ambiguity) '}' (rule end)
	 *     labels+=ElkLabel ','? ']' (ambiguity) '}' (rule end)
	 *     properties+=Property ','? '}' (ambiguity) '}' (rule end)
	 *     sources+=[ElkConnectableShape|STRING] ','? ']' (ambiguity) '}' (rule end)
	 *     targets+=[ElkConnectableShape|STRING] ','? ']' (ambiguity) '}' (rule end)
	 */
	protected void emit_ElkEdge_CommaKeyword_3_q(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
	/**
	 * Ambiguous syntax:
	 *     ','?
	 *
	 * This ambiguous syntax occurs at:
	 *     labels+=ElkLabel (ambiguity) ']' ',' KeyChildren ':' '[' children+=ElkNode
	 *     labels+=ElkLabel (ambiguity) ']' ',' KeyEdges ':' '[' containedEdges+=ElkEdge
	 *     labels+=ElkLabel (ambiguity) ']' ',' KeyHeight ':' height=Number
	 *     labels+=ElkLabel (ambiguity) ']' ',' KeyId ':' identifier=STRING
	 *     labels+=ElkLabel (ambiguity) ']' ',' KeyLabels ':' '[' labels+=ElkLabel
	 *     labels+=ElkLabel (ambiguity) ']' ',' KeyLayoutOptions ':' '{' properties+=Property
	 *     labels+=ElkLabel (ambiguity) ']' ',' KeyPorts ':' '[' ports+=ElkPort
	 *     labels+=ElkLabel (ambiguity) ']' ',' KeySources ':' '[' sources+=[ElkConnectableShape|STRING]
	 *     labels+=ElkLabel (ambiguity) ']' ',' KeyTargets ':' '[' targets+=[ElkConnectableShape|STRING]
	 *     labels+=ElkLabel (ambiguity) ']' ',' KeyText ':' text=STRING
	 *     labels+=ElkLabel (ambiguity) ']' ',' KeyWidth ':' width=Number
	 *     labels+=ElkLabel (ambiguity) ']' ',' KeyX ':' x=Number
	 *     labels+=ElkLabel (ambiguity) ']' ',' KeyY ':' y=Number
	 *     labels+=ElkLabel (ambiguity) ']' ','? '}' (rule end)
	 */
	protected void emit_ElkGraphElementLabels_CommaKeyword_2_q(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
	/**
	 * Ambiguous syntax:
	 *     (
	  *         ','? 
	  *         ']' 
	  *         ',' 
	  *         KeyLabels 
	  *         ':' 
	  *         '['
	  *     )+
	 *
	 * This ambiguous syntax occurs at:
	 *     labels+=ElkLabel (ambiguity) labels+=ElkLabel
	 */
	protected void emit_ElkGraphElementLabels_ElkLabel_LabelElement___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_0_KeyLabelsParserRuleCall_3_0_ColonKeyword_3_1_LeftSquareBracketKeyword_0__p(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
	/**
	 * Ambiguous syntax:
	 *     (
	  *         ','? 
	  *         ']' 
	  *         ',' 
	  *         KeyLabels 
	  *         ':' 
	  *         '['
	  *     )+
	 *
	 * This ambiguous syntax occurs at:
	 *     labels+=ElkLabel (ambiguity) labels+=ElkLabel
	 */
	protected void emit_ElkGraphElementLabels_ElkNode_NodeElement___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_1_0_KeyLabelsParserRuleCall_4_0_ColonKeyword_4_1_LeftSquareBracketKeyword_0__p(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
	/**
	 * Ambiguous syntax:
	 *     (
	  *         ','? 
	  *         ']' 
	  *         ',' 
	  *         KeyLabels 
	  *         ':' 
	  *         '['
	  *     )+
	 *
	 * This ambiguous syntax occurs at:
	 *     labels+=ElkLabel (ambiguity) labels+=ElkLabel
	 */
	protected void emit_ElkGraphElementLabels_ElkPort_PortElement___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_0_KeyLabelsParserRuleCall_2_0_ColonKeyword_2_1_LeftSquareBracketKeyword_0__p(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
	/**
	 * Ambiguous syntax:
	 *     ','?
	 *
	 * This ambiguous syntax occurs at:
	 *     properties+=Property (ambiguity) '}' ',' KeyChildren ':' '[' children+=ElkNode
	 *     properties+=Property (ambiguity) '}' ',' KeyEdges ':' '[' containedEdges+=ElkEdge
	 *     properties+=Property (ambiguity) '}' ',' KeyHeight ':' height=Number
	 *     properties+=Property (ambiguity) '}' ',' KeyId ':' identifier=STRING
	 *     properties+=Property (ambiguity) '}' ',' KeyLabels ':' '[' labels+=ElkLabel
	 *     properties+=Property (ambiguity) '}' ',' KeyLayoutOptions ':' '{' properties+=Property
	 *     properties+=Property (ambiguity) '}' ',' KeyPorts ':' '[' ports+=ElkPort
	 *     properties+=Property (ambiguity) '}' ',' KeySources ':' '[' sources+=[ElkConnectableShape|STRING]
	 *     properties+=Property (ambiguity) '}' ',' KeyTargets ':' '[' targets+=[ElkConnectableShape|STRING]
	 *     properties+=Property (ambiguity) '}' ',' KeyText ':' text=STRING
	 *     properties+=Property (ambiguity) '}' ',' KeyWidth ':' width=Number
	 *     properties+=Property (ambiguity) '}' ',' KeyX ':' x=Number
	 *     properties+=Property (ambiguity) '}' ',' KeyY ':' y=Number
	 *     properties+=Property (ambiguity) '}' ','? '}' (rule end)
	 */
	protected void emit_ElkGraphElementProperties_CommaKeyword_2_q(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
	/**
	 * Ambiguous syntax:
	 *     (
	  *         ','? 
	  *         '}' 
	  *         ',' 
	  *         KeyLayoutOptions 
	  *         ':' 
	  *         '{'
	  *     )+
	 *
	 * This ambiguous syntax occurs at:
	 *     properties+=Property (ambiguity) properties+=Property
	 */
	protected void emit_ElkGraphElementProperties_ElkLabel_LabelElement___CommaKeyword_2_q_RightCurlyBracketKeyword_3_CommaKeyword_2_0_KeyLayoutOptionsParserRuleCall_4_0_ColonKeyword_4_1_LeftCurlyBracketKeyword_0__p(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
	/**
	 * Ambiguous syntax:
	 *     (
	  *         ','? 
	  *         '}' 
	  *         ',' 
	  *         KeyLayoutOptions 
	  *         ':' 
	  *         '{'
	  *     )+
	 *
	 * This ambiguous syntax occurs at:
	 *     properties+=Property (ambiguity) properties+=Property
	 */
	protected void emit_ElkGraphElementProperties_ElkNode_NodeElement___CommaKeyword_2_q_RightCurlyBracketKeyword_3_CommaKeyword_2_1_0_KeyLayoutOptionsParserRuleCall_6_0_ColonKeyword_6_1_LeftCurlyBracketKeyword_0__p(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
	/**
	 * Ambiguous syntax:
	 *     (
	  *         ','? 
	  *         '}' 
	  *         ',' 
	  *         KeyLayoutOptions 
	  *         ':' 
	  *         '{'
	  *     )+
	 *
	 * This ambiguous syntax occurs at:
	 *     properties+=Property (ambiguity) properties+=Property
	 */
	protected void emit_ElkGraphElementProperties_ElkPort_PortElement___CommaKeyword_2_q_RightCurlyBracketKeyword_3_CommaKeyword_2_0_KeyLayoutOptionsParserRuleCall_3_0_ColonKeyword_3_1_LeftCurlyBracketKeyword_0__p(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
	/**
	 * Ambiguous syntax:
	 *     ','?
	 *
	 * This ambiguous syntax occurs at:
	 *     (rule start) '{' JsonMember (ambiguity) '}' (rule start)
	 *     height=Number (ambiguity) '}' (rule end)
	 *     identifier=STRING (ambiguity) '}' (rule end)
	 *     labels+=ElkLabel ','? ']' (ambiguity) '}' (rule end)
	 *     properties+=Property ','? '}' (ambiguity) '}' (rule end)
	 *     text=STRING (ambiguity) '}' (rule end)
	 *     width=Number (ambiguity) '}' (rule end)
	 *     x=Number (ambiguity) '}' (rule end)
	 *     y=Number (ambiguity) '}' (rule end)
	 */
	protected void emit_ElkLabel_CommaKeyword_3_q(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
	/**
	 * Ambiguous syntax:
	 *     (JsonMember ',')*
	 *
	 * This ambiguous syntax occurs at:
	 *     text=STRING ',' (ambiguity) KeyText ':' text=STRING
	 */
	protected void emit_ElkLabel_LabelElement___JsonMemberParserRuleCall_5_CommaKeyword_2_0__a(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
	/**
	 * Ambiguous syntax:
	 *     ','?
	 *
	 * This ambiguous syntax occurs at:
	 *     children+=ElkNode (ambiguity) ']' ',' KeyChildren ':' '[' children+=ElkNode
	 *     children+=ElkNode (ambiguity) ']' ',' KeyEdges ':' '[' containedEdges+=ElkEdge
	 *     children+=ElkNode (ambiguity) ']' ',' KeyHeight ':' height=Number
	 *     children+=ElkNode (ambiguity) ']' ',' KeyId ':' identifier=STRING
	 *     children+=ElkNode (ambiguity) ']' ',' KeyLabels ':' '[' labels+=ElkLabel
	 *     children+=ElkNode (ambiguity) ']' ',' KeyLayoutOptions ':' '{' properties+=Property
	 *     children+=ElkNode (ambiguity) ']' ',' KeyPorts ':' '[' ports+=ElkPort
	 *     children+=ElkNode (ambiguity) ']' ',' KeyWidth ':' width=Number
	 *     children+=ElkNode (ambiguity) ']' ',' KeyX ':' x=Number
	 *     children+=ElkNode (ambiguity) ']' ',' KeyY ':' y=Number
	 *     children+=ElkNode (ambiguity) ']' ','? '}' (rule end)
	 */
	protected void emit_ElkNodeChildren_CommaKeyword_2_q(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
	/**
	 * Ambiguous syntax:
	 *     ','?
	 *
	 * This ambiguous syntax occurs at:
	 *     containedEdges+=ElkEdge (ambiguity) ']' ',' KeyChildren ':' '[' children+=ElkNode
	 *     containedEdges+=ElkEdge (ambiguity) ']' ',' KeyEdges ':' '[' containedEdges+=ElkEdge
	 *     containedEdges+=ElkEdge (ambiguity) ']' ',' KeyHeight ':' height=Number
	 *     containedEdges+=ElkEdge (ambiguity) ']' ',' KeyId ':' identifier=STRING
	 *     containedEdges+=ElkEdge (ambiguity) ']' ',' KeyLabels ':' '[' labels+=ElkLabel
	 *     containedEdges+=ElkEdge (ambiguity) ']' ',' KeyLayoutOptions ':' '{' properties+=Property
	 *     containedEdges+=ElkEdge (ambiguity) ']' ',' KeyPorts ':' '[' ports+=ElkPort
	 *     containedEdges+=ElkEdge (ambiguity) ']' ',' KeyWidth ':' width=Number
	 *     containedEdges+=ElkEdge (ambiguity) ']' ',' KeyX ':' x=Number
	 *     containedEdges+=ElkEdge (ambiguity) ']' ',' KeyY ':' y=Number
	 *     containedEdges+=ElkEdge (ambiguity) ']' ','? '}' (rule end)
	 */
	protected void emit_ElkNodeEdges_CommaKeyword_2_q(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
	/**
	 * Ambiguous syntax:
	 *     ','?
	 *
	 * This ambiguous syntax occurs at:
	 *     ports+=ElkPort (ambiguity) ']' ',' KeyChildren ':' '[' children+=ElkNode
	 *     ports+=ElkPort (ambiguity) ']' ',' KeyEdges ':' '[' containedEdges+=ElkEdge
	 *     ports+=ElkPort (ambiguity) ']' ',' KeyHeight ':' height=Number
	 *     ports+=ElkPort (ambiguity) ']' ',' KeyId ':' identifier=STRING
	 *     ports+=ElkPort (ambiguity) ']' ',' KeyLabels ':' '[' labels+=ElkLabel
	 *     ports+=ElkPort (ambiguity) ']' ',' KeyLayoutOptions ':' '{' properties+=Property
	 *     ports+=ElkPort (ambiguity) ']' ',' KeyPorts ':' '[' ports+=ElkPort
	 *     ports+=ElkPort (ambiguity) ']' ',' KeyWidth ':' width=Number
	 *     ports+=ElkPort (ambiguity) ']' ',' KeyX ':' x=Number
	 *     ports+=ElkPort (ambiguity) ']' ',' KeyY ':' y=Number
	 *     ports+=ElkPort (ambiguity) ']' ','? '}' (rule end)
	 */
	protected void emit_ElkNodePorts_CommaKeyword_2_q(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
	/**
	 * Ambiguous syntax:
	 *     ','?
	 *
	 * This ambiguous syntax occurs at:
	 *     (rule start) '{' (ambiguity) '}' (rule start)
	 *     children+=ElkNode ','? ']' (ambiguity) '}' (rule end)
	 *     containedEdges+=ElkEdge ','? ']' (ambiguity) '}' (rule end)
	 *     height=Number (ambiguity) '}' (rule end)
	 *     identifier=STRING (ambiguity) '}' (rule end)
	 *     labels+=ElkLabel ','? ']' (ambiguity) '}' (rule end)
	 *     ports+=ElkPort ','? ']' (ambiguity) '}' (rule end)
	 *     properties+=Property ','? '}' (ambiguity) '}' (rule end)
	 *     width=Number (ambiguity) '}' (rule end)
	 *     x=Number (ambiguity) '}' (rule end)
	 *     y=Number (ambiguity) '}' (rule end)
	 */
	protected void emit_ElkNode_CommaKeyword_3_q(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
	/**
	 * Ambiguous syntax:
	 *     (
	  *         ','? 
	  *         ']' 
	  *         ',' 
	  *         KeyChildren 
	  *         ':' 
	  *         '['
	  *     )+
	 *
	 * This ambiguous syntax occurs at:
	 *     children+=ElkNode (ambiguity) children+=ElkNode
	 */
	protected void emit_ElkNode_ElkNodeChildren_NodeElement___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_1_0_KeyChildrenParserRuleCall_2_0_ColonKeyword_2_1_LeftSquareBracketKeyword_0__p(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
	/**
	 * Ambiguous syntax:
	 *     (
	  *         ','? 
	  *         ']' 
	  *         ',' 
	  *         KeyEdges 
	  *         ':' 
	  *         '['
	  *     )+
	 *
	 * This ambiguous syntax occurs at:
	 *     containedEdges+=ElkEdge (ambiguity) containedEdges+=ElkEdge
	 */
	protected void emit_ElkNode_ElkNodeEdges_NodeElement___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_1_0_KeyEdgesParserRuleCall_5_0_ColonKeyword_5_1_LeftSquareBracketKeyword_0__p(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
	/**
	 * Ambiguous syntax:
	 *     (
	  *         ','? 
	  *         ']' 
	  *         ',' 
	  *         KeyPorts 
	  *         ':' 
	  *         '['
	  *     )+
	 *
	 * This ambiguous syntax occurs at:
	 *     ports+=ElkPort (ambiguity) ports+=ElkPort
	 */
	protected void emit_ElkNode_ElkNodePorts_NodeElement___CommaKeyword_2_q_RightSquareBracketKeyword_3_CommaKeyword_2_1_0_KeyPortsParserRuleCall_3_0_ColonKeyword_3_1_LeftSquareBracketKeyword_0__p(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
	/**
	 * Ambiguous syntax:
	 *     ','?
	 *
	 * This ambiguous syntax occurs at:
	 *     (rule start) '{' JsonMember (ambiguity) '}' (rule start)
	 *     height=Number (ambiguity) '}' (rule end)
	 *     identifier=STRING (ambiguity) '}' (rule end)
	 *     labels+=ElkLabel ','? ']' (ambiguity) '}' (rule end)
	 *     properties+=Property ','? '}' (ambiguity) '}' (rule end)
	 *     width=Number (ambiguity) '}' (rule end)
	 *     x=Number (ambiguity) '}' (rule end)
	 *     y=Number (ambiguity) '}' (rule end)
	 */
	protected void emit_ElkPort_CommaKeyword_3_q(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
}
