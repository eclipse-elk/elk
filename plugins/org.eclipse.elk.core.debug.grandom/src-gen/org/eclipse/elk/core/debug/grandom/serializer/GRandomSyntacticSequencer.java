/*******************************************************************************
 * Copyright (c) 2018 Kiel University and others.
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
package org.eclipse.elk.core.debug.grandom.serializer;

import com.google.inject.Inject;
import java.util.List;
import org.eclipse.elk.core.debug.grandom.services.GRandomGrammarAccess;
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
public class GRandomSyntacticSequencer extends AbstractSyntacticSequencer {

	protected GRandomGrammarAccess grammarAccess;
	protected AbstractElementAlias match_Configuration___LeftCurlyBracketKeyword_3_0_RightCurlyBracketKeyword_3_2__q;
	protected AbstractElementAlias match_Edges___LeftCurlyBracketKeyword_1_0_RightCurlyBracketKeyword_1_2__q;
	protected AbstractElementAlias match_Hierarchy___LeftCurlyBracketKeyword_2_0_RightCurlyBracketKeyword_2_2__q;
	protected AbstractElementAlias match_Nodes___LeftCurlyBracketKeyword_4_0_RightCurlyBracketKeyword_4_2__q;
	protected AbstractElementAlias match_Ports___LeftCurlyBracketKeyword_2_0_RightCurlyBracketKeyword_2_2__q;
	protected AbstractElementAlias match_Size___LeftCurlyBracketKeyword_1_1_0_RightCurlyBracketKeyword_1_1_2__q;
	
	@Inject
	protected void init(IGrammarAccess access) {
		grammarAccess = (GRandomGrammarAccess) access;
		match_Configuration___LeftCurlyBracketKeyword_3_0_RightCurlyBracketKeyword_3_2__q = new GroupAlias(false, true, new TokenAlias(false, false, grammarAccess.getConfigurationAccess().getLeftCurlyBracketKeyword_3_0()), new TokenAlias(false, false, grammarAccess.getConfigurationAccess().getRightCurlyBracketKeyword_3_2()));
		match_Edges___LeftCurlyBracketKeyword_1_0_RightCurlyBracketKeyword_1_2__q = new GroupAlias(false, true, new TokenAlias(false, false, grammarAccess.getEdgesAccess().getLeftCurlyBracketKeyword_1_0()), new TokenAlias(false, false, grammarAccess.getEdgesAccess().getRightCurlyBracketKeyword_1_2()));
		match_Hierarchy___LeftCurlyBracketKeyword_2_0_RightCurlyBracketKeyword_2_2__q = new GroupAlias(false, true, new TokenAlias(false, false, grammarAccess.getHierarchyAccess().getLeftCurlyBracketKeyword_2_0()), new TokenAlias(false, false, grammarAccess.getHierarchyAccess().getRightCurlyBracketKeyword_2_2()));
		match_Nodes___LeftCurlyBracketKeyword_4_0_RightCurlyBracketKeyword_4_2__q = new GroupAlias(false, true, new TokenAlias(false, false, grammarAccess.getNodesAccess().getLeftCurlyBracketKeyword_4_0()), new TokenAlias(false, false, grammarAccess.getNodesAccess().getRightCurlyBracketKeyword_4_2()));
		match_Ports___LeftCurlyBracketKeyword_2_0_RightCurlyBracketKeyword_2_2__q = new GroupAlias(false, true, new TokenAlias(false, false, grammarAccess.getPortsAccess().getLeftCurlyBracketKeyword_2_0()), new TokenAlias(false, false, grammarAccess.getPortsAccess().getRightCurlyBracketKeyword_2_2()));
		match_Size___LeftCurlyBracketKeyword_1_1_0_RightCurlyBracketKeyword_1_1_2__q = new GroupAlias(false, true, new TokenAlias(false, false, grammarAccess.getSizeAccess().getLeftCurlyBracketKeyword_1_1_0()), new TokenAlias(false, false, grammarAccess.getSizeAccess().getRightCurlyBracketKeyword_1_1_2()));
	}
	
	@Override
	protected String getUnassignedRuleCallToken(EObject semanticObject, RuleCall ruleCall, INode node) {
		if (ruleCall.getRule() == grammarAccess.getLabelsRule())
			return getLabelsToken(semanticObject, ruleCall, node);
		else if (ruleCall.getRule() == grammarAccess.getPmRule())
			return getPmToken(semanticObject, ruleCall, node);
		return "";
	}
	
	/**
	 * Labels:
	 *     'labels'
	 * ;
	 */
	protected String getLabelsToken(EObject semanticObject, RuleCall ruleCall, INode node) {
		if (node != null)
			return getTokenText(node);
		return "labels";
	}
	
	/**
	 * Pm:
	 *     '+/-'
	 * ;
	 */
	protected String getPmToken(EObject semanticObject, RuleCall ruleCall, INode node) {
		if (node != null)
			return getTokenText(node);
		return "+/-";
	}
	
	@Override
	protected void emitUnassignedTokens(EObject semanticObject, ISynTransition transition, INode fromNode, INode toNode) {
		if (transition.getAmbiguousSyntaxes().isEmpty()) return;
		List<INode> transitionNodes = collectNodes(fromNode, toNode);
		for (AbstractElementAlias syntax : transition.getAmbiguousSyntaxes()) {
			List<INode> syntaxNodes = getNodesFor(transitionNodes, syntax);
			if (match_Configuration___LeftCurlyBracketKeyword_3_0_RightCurlyBracketKeyword_3_2__q.equals(syntax))
				emit_Configuration___LeftCurlyBracketKeyword_3_0_RightCurlyBracketKeyword_3_2__q(semanticObject, getLastNavigableState(), syntaxNodes);
			else if (match_Edges___LeftCurlyBracketKeyword_1_0_RightCurlyBracketKeyword_1_2__q.equals(syntax))
				emit_Edges___LeftCurlyBracketKeyword_1_0_RightCurlyBracketKeyword_1_2__q(semanticObject, getLastNavigableState(), syntaxNodes);
			else if (match_Hierarchy___LeftCurlyBracketKeyword_2_0_RightCurlyBracketKeyword_2_2__q.equals(syntax))
				emit_Hierarchy___LeftCurlyBracketKeyword_2_0_RightCurlyBracketKeyword_2_2__q(semanticObject, getLastNavigableState(), syntaxNodes);
			else if (match_Nodes___LeftCurlyBracketKeyword_4_0_RightCurlyBracketKeyword_4_2__q.equals(syntax))
				emit_Nodes___LeftCurlyBracketKeyword_4_0_RightCurlyBracketKeyword_4_2__q(semanticObject, getLastNavigableState(), syntaxNodes);
			else if (match_Ports___LeftCurlyBracketKeyword_2_0_RightCurlyBracketKeyword_2_2__q.equals(syntax))
				emit_Ports___LeftCurlyBracketKeyword_2_0_RightCurlyBracketKeyword_2_2__q(semanticObject, getLastNavigableState(), syntaxNodes);
			else if (match_Size___LeftCurlyBracketKeyword_1_1_0_RightCurlyBracketKeyword_1_1_2__q.equals(syntax))
				emit_Size___LeftCurlyBracketKeyword_1_1_0_RightCurlyBracketKeyword_1_1_2__q(semanticObject, getLastNavigableState(), syntaxNodes);
			else acceptNodes(getLastNavigableState(), syntaxNodes);
		}
	}

	/**
	 * Ambiguous syntax:
	 *     ('{' '}')?
	 *
	 * This ambiguous syntax occurs at:
	 *     form=Form (ambiguity) (rule end)
	 */
	protected void emit_Configuration___LeftCurlyBracketKeyword_3_0_RightCurlyBracketKeyword_3_2__q(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
	/**
	 * Ambiguous syntax:
	 *     ('{' '}')?
	 *
	 * This ambiguous syntax occurs at:
	 *     nEdges=DoubleQuantity (ambiguity) (rule end)
	 */
	protected void emit_Edges___LeftCurlyBracketKeyword_1_0_RightCurlyBracketKeyword_1_2__q(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
	/**
	 * Ambiguous syntax:
	 *     ('{' '}')?
	 *
	 * This ambiguous syntax occurs at:
	 *     (rule start) 'hierarchy' (ambiguity) (rule start)
	 */
	protected void emit_Hierarchy___LeftCurlyBracketKeyword_2_0_RightCurlyBracketKeyword_2_2__q(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
	/**
	 * Ambiguous syntax:
	 *     ('{' '}')?
	 *
	 * This ambiguous syntax occurs at:
	 *     nNodes=DoubleQuantity (ambiguity) (rule end)
	 */
	protected void emit_Nodes___LeftCurlyBracketKeyword_4_0_RightCurlyBracketKeyword_4_2__q(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
	/**
	 * Ambiguous syntax:
	 *     ('{' '}')?
	 *
	 * This ambiguous syntax occurs at:
	 *     (rule start) 'ports' (ambiguity) (rule start)
	 */
	protected void emit_Ports___LeftCurlyBracketKeyword_2_0_RightCurlyBracketKeyword_2_2__q(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
	/**
	 * Ambiguous syntax:
	 *     ('{' '}')?
	 *
	 * This ambiguous syntax occurs at:
	 *     (rule start) 'size' (ambiguity) (rule start)
	 */
	protected void emit_Size___LeftCurlyBracketKeyword_1_1_0_RightCurlyBracketKeyword_1_1_2__q(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
}
