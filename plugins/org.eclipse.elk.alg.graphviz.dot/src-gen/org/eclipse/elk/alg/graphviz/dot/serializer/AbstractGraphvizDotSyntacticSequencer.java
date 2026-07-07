/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.graphviz.dot.serializer;

import com.google.inject.Inject;
import java.util.List;
import org.eclipse.elk.alg.graphviz.dot.services.GraphvizDotGrammarAccess;
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
public abstract class AbstractGraphvizDotSyntacticSequencer extends AbstractSyntacticSequencer {

	protected GraphvizDotGrammarAccess grammarAccess;
	protected AbstractElementAlias match_AttributeStatement_CommaKeyword_2_1_0_q;
	protected AbstractElementAlias match_EdgeStatement_CommaKeyword_2_1_1_0_q;
	protected AbstractElementAlias match_EdgeStatement___LeftSquareBracketKeyword_2_0_RightSquareBracketKeyword_2_2__q;
	protected AbstractElementAlias match_NodeStatement_CommaKeyword_1_1_1_0_q;
	protected AbstractElementAlias match_NodeStatement___LeftSquareBracketKeyword_1_0_RightSquareBracketKeyword_1_2__q;
	protected AbstractElementAlias match_Statement_SemicolonKeyword_1_q;
	protected AbstractElementAlias match_Subgraph_SubgraphKeyword_1_0_q;
	
	@Inject
	protected void init(IGrammarAccess access) {
		grammarAccess = (GraphvizDotGrammarAccess) access;
		match_AttributeStatement_CommaKeyword_2_1_0_q = new TokenAlias(false, true, grammarAccess.getAttributeStatementAccess().getCommaKeyword_2_1_0());
		match_EdgeStatement_CommaKeyword_2_1_1_0_q = new TokenAlias(false, true, grammarAccess.getEdgeStatementAccess().getCommaKeyword_2_1_1_0());
		match_EdgeStatement___LeftSquareBracketKeyword_2_0_RightSquareBracketKeyword_2_2__q = new GroupAlias(false, true, new TokenAlias(false, false, grammarAccess.getEdgeStatementAccess().getLeftSquareBracketKeyword_2_0()), new TokenAlias(false, false, grammarAccess.getEdgeStatementAccess().getRightSquareBracketKeyword_2_2()));
		match_NodeStatement_CommaKeyword_1_1_1_0_q = new TokenAlias(false, true, grammarAccess.getNodeStatementAccess().getCommaKeyword_1_1_1_0());
		match_NodeStatement___LeftSquareBracketKeyword_1_0_RightSquareBracketKeyword_1_2__q = new GroupAlias(false, true, new TokenAlias(false, false, grammarAccess.getNodeStatementAccess().getLeftSquareBracketKeyword_1_0()), new TokenAlias(false, false, grammarAccess.getNodeStatementAccess().getRightSquareBracketKeyword_1_2()));
		match_Statement_SemicolonKeyword_1_q = new TokenAlias(false, true, grammarAccess.getStatementAccess().getSemicolonKeyword_1());
		match_Subgraph_SubgraphKeyword_1_0_q = new TokenAlias(false, true, grammarAccess.getSubgraphAccess().getSubgraphKeyword_1_0());
	}
	
	@Override
	protected String getUnassignedRuleCallToken(EObject semanticObject, RuleCall ruleCall, INode node) {
		return "";
	}
	
	
	@Override
	protected void emitUnassignedTokens(EObject semanticObject, ISynTransition transition, INode fromNode, INode toNode) {
		if (transition.getAmbiguousSyntaxes().isEmpty()) return;
		List<INode> transitionNodes = collectNodes(fromNode, toNode);
		for (AbstractElementAlias syntax : transition.getAmbiguousSyntaxes()) {
			List<INode> syntaxNodes = getNodesFor(transitionNodes, syntax);
			if (match_AttributeStatement_CommaKeyword_2_1_0_q.equals(syntax))
				emit_AttributeStatement_CommaKeyword_2_1_0_q(semanticObject, getLastNavigableState(), syntaxNodes);
			else if (match_EdgeStatement_CommaKeyword_2_1_1_0_q.equals(syntax))
				emit_EdgeStatement_CommaKeyword_2_1_1_0_q(semanticObject, getLastNavigableState(), syntaxNodes);
			else if (match_EdgeStatement___LeftSquareBracketKeyword_2_0_RightSquareBracketKeyword_2_2__q.equals(syntax))
				emit_EdgeStatement___LeftSquareBracketKeyword_2_0_RightSquareBracketKeyword_2_2__q(semanticObject, getLastNavigableState(), syntaxNodes);
			else if (match_NodeStatement_CommaKeyword_1_1_1_0_q.equals(syntax))
				emit_NodeStatement_CommaKeyword_1_1_1_0_q(semanticObject, getLastNavigableState(), syntaxNodes);
			else if (match_NodeStatement___LeftSquareBracketKeyword_1_0_RightSquareBracketKeyword_1_2__q.equals(syntax))
				emit_NodeStatement___LeftSquareBracketKeyword_1_0_RightSquareBracketKeyword_1_2__q(semanticObject, getLastNavigableState(), syntaxNodes);
			else if (match_Statement_SemicolonKeyword_1_q.equals(syntax))
				emit_Statement_SemicolonKeyword_1_q(semanticObject, getLastNavigableState(), syntaxNodes);
			else if (match_Subgraph_SubgraphKeyword_1_0_q.equals(syntax))
				emit_Subgraph_SubgraphKeyword_1_0_q(semanticObject, getLastNavigableState(), syntaxNodes);
			else acceptNodes(getLastNavigableState(), syntaxNodes);
		}
	}

	/**
	 * Ambiguous syntax:
	 *     ','?
	 *
	 * This ambiguous syntax occurs at:
	 *     attributes+=ListAttribute (ambiguity) attributes+=ListAttribute
	 */
	protected void emit_AttributeStatement_CommaKeyword_2_1_0_q(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
	/**
	 * Ambiguous syntax:
	 *     ','?
	 *
	 * This ambiguous syntax occurs at:
	 *     attributes+=ListAttribute (ambiguity) attributes+=ListAttribute
	 */
	protected void emit_EdgeStatement_CommaKeyword_2_1_1_0_q(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
	/**
	 * Ambiguous syntax:
	 *     ('[' ']')?
	 *
	 * This ambiguous syntax occurs at:
	 *     edgeTargets+=EdgeTarget (ambiguity) ';'? (rule end)
	 *     edgeTargets+=EdgeTarget (ambiguity) (rule end)
	 */
	protected void emit_EdgeStatement___LeftSquareBracketKeyword_2_0_RightSquareBracketKeyword_2_2__q(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
	/**
	 * Ambiguous syntax:
	 *     ','?
	 *
	 * This ambiguous syntax occurs at:
	 *     attributes+=ListAttribute (ambiguity) attributes+=ListAttribute
	 */
	protected void emit_NodeStatement_CommaKeyword_1_1_1_0_q(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
	/**
	 * Ambiguous syntax:
	 *     ('[' ']')?
	 *
	 * This ambiguous syntax occurs at:
	 *     node=Node (ambiguity) ';'? (rule end)
	 *     node=Node (ambiguity) (rule end)
	 */
	protected void emit_NodeStatement___LeftSquareBracketKeyword_1_0_RightSquareBracketKeyword_1_2__q(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
	/**
	 * Ambiguous syntax:
	 *     ';'?
	 *
	 * This ambiguous syntax occurs at:
	 *     (rule start) 'subgraph'? '{' '}' (ambiguity) (rule start)
	 *     attributes+=ListAttribute ']' (ambiguity) (rule end)
	 *     edgeTargets+=EdgeTarget ('[' ']')? (ambiguity) (rule end)
	 *     name=ID '{' '}' (ambiguity) (rule end)
	 *     node=Node ('[' ']')? (ambiguity) (rule end)
	 *     statements+=Statement '}' (ambiguity) (rule end)
	 *     type=AttributeType '[' ']' (ambiguity) (rule end)
	 *     value=DotID (ambiguity) (rule end)
	 */
	protected void emit_Statement_SemicolonKeyword_1_q(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
	/**
	 * Ambiguous syntax:
	 *     'subgraph'?
	 *
	 * This ambiguous syntax occurs at:
	 *     (rule start) (ambiguity) '{' '}' ';'? (rule start)
	 *     (rule start) (ambiguity) '{' '}' (rule start)
	 *     (rule start) (ambiguity) '{' statements+=Statement
	 */
	protected void emit_Subgraph_SubgraphKeyword_1_0_q(EObject semanticObject, ISynNavigable transition, List<INode> nodes) {
		acceptNodes(transition, nodes);
	}
	
}
