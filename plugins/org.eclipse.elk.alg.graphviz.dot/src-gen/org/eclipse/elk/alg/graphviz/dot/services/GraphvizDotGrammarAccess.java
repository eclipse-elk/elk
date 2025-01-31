/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.graphviz.dot.services;

import com.google.inject.Inject;
import com.google.inject.Singleton;
import java.util.List;
import org.eclipse.xtext.Action;
import org.eclipse.xtext.Alternatives;
import org.eclipse.xtext.Assignment;
import org.eclipse.xtext.EnumLiteralDeclaration;
import org.eclipse.xtext.EnumRule;
import org.eclipse.xtext.Grammar;
import org.eclipse.xtext.GrammarUtil;
import org.eclipse.xtext.Group;
import org.eclipse.xtext.Keyword;
import org.eclipse.xtext.ParserRule;
import org.eclipse.xtext.RuleCall;
import org.eclipse.xtext.TerminalRule;
import org.eclipse.xtext.common.services.TerminalsGrammarAccess;
import org.eclipse.xtext.service.AbstractElementFinder.AbstractGrammarElementFinder;
import org.eclipse.xtext.service.GrammarProvider;

@Singleton
public class GraphvizDotGrammarAccess extends AbstractGrammarElementFinder {
	
	public class GraphvizModelElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.alg.graphviz.dot.GraphvizDot.GraphvizModel");
		private final Assignment cGraphsAssignment = (Assignment)rule.eContents().get(1);
		private final RuleCall cGraphsGraphParserRuleCall_0 = (RuleCall)cGraphsAssignment.eContents().get(0);
		
		//GraphvizModel:
		//	graphs+=Graph*;
		@Override public ParserRule getRule() { return rule; }
		
		//graphs+=Graph*
		public Assignment getGraphsAssignment() { return cGraphsAssignment; }
		
		//Graph
		public RuleCall getGraphsGraphParserRuleCall_0() { return cGraphsGraphParserRuleCall_0; }
	}
	public class GraphElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.alg.graphviz.dot.GraphvizDot.Graph");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Assignment cStrictAssignment_0 = (Assignment)cGroup.eContents().get(0);
		private final Keyword cStrictStrictKeyword_0_0 = (Keyword)cStrictAssignment_0.eContents().get(0);
		private final Assignment cTypeAssignment_1 = (Assignment)cGroup.eContents().get(1);
		private final RuleCall cTypeGraphTypeEnumRuleCall_1_0 = (RuleCall)cTypeAssignment_1.eContents().get(0);
		private final Assignment cNameAssignment_2 = (Assignment)cGroup.eContents().get(2);
		private final RuleCall cNameDotIDParserRuleCall_2_0 = (RuleCall)cNameAssignment_2.eContents().get(0);
		private final Keyword cLeftCurlyBracketKeyword_3 = (Keyword)cGroup.eContents().get(3);
		private final Assignment cStatementsAssignment_4 = (Assignment)cGroup.eContents().get(4);
		private final RuleCall cStatementsStatementParserRuleCall_4_0 = (RuleCall)cStatementsAssignment_4.eContents().get(0);
		private final Keyword cRightCurlyBracketKeyword_5 = (Keyword)cGroup.eContents().get(5);
		
		///*~~~~Original definition~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// *  graph:      [ 'strict' ] ('graph' | 'digraph') [ ID ] '{' stmt_list '}'
		// */ Graph:
		//	strict?="strict"? type=GraphType name=DotID?
		//	"{" statements+=Statement* "}";
		@Override public ParserRule getRule() { return rule; }
		
		//strict?="strict"? type=GraphType name=DotID? "{" statements+=Statement* "}"
		public Group getGroup() { return cGroup; }
		
		//strict?="strict"?
		public Assignment getStrictAssignment_0() { return cStrictAssignment_0; }
		
		//"strict"
		public Keyword getStrictStrictKeyword_0_0() { return cStrictStrictKeyword_0_0; }
		
		//type=GraphType
		public Assignment getTypeAssignment_1() { return cTypeAssignment_1; }
		
		//GraphType
		public RuleCall getTypeGraphTypeEnumRuleCall_1_0() { return cTypeGraphTypeEnumRuleCall_1_0; }
		
		//name=DotID?
		public Assignment getNameAssignment_2() { return cNameAssignment_2; }
		
		//DotID
		public RuleCall getNameDotIDParserRuleCall_2_0() { return cNameDotIDParserRuleCall_2_0; }
		
		//"{"
		public Keyword getLeftCurlyBracketKeyword_3() { return cLeftCurlyBracketKeyword_3; }
		
		//statements+=Statement*
		public Assignment getStatementsAssignment_4() { return cStatementsAssignment_4; }
		
		//Statement
		public RuleCall getStatementsStatementParserRuleCall_4_0() { return cStatementsStatementParserRuleCall_4_0; }
		
		//"}"
		public Keyword getRightCurlyBracketKeyword_5() { return cRightCurlyBracketKeyword_5; }
	}
	public class StatementElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.alg.graphviz.dot.GraphvizDot.Statement");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Alternatives cAlternatives_0 = (Alternatives)cGroup.eContents().get(0);
		private final RuleCall cNodeStatementParserRuleCall_0_0 = (RuleCall)cAlternatives_0.eContents().get(0);
		private final RuleCall cEdgeStatementParserRuleCall_0_1 = (RuleCall)cAlternatives_0.eContents().get(1);
		private final RuleCall cAttributeStatementParserRuleCall_0_2 = (RuleCall)cAlternatives_0.eContents().get(2);
		private final RuleCall cAttributeParserRuleCall_0_3 = (RuleCall)cAlternatives_0.eContents().get(3);
		private final RuleCall cSubgraphParserRuleCall_0_4 = (RuleCall)cAlternatives_0.eContents().get(4);
		private final Keyword cSemicolonKeyword_1 = (Keyword)cGroup.eContents().get(1);
		
		///*~~~~Original definition~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// *  stmt_list:  [ stmt [ ';' ] [ stmt_list ] ]
		// *  stmt:       node_stmt | edge_stmt | attr_stmt | ID '=' ID | subgraph
		// */ Statement:
		//	(NodeStatement | EdgeStatement | AttributeStatement | Attribute
		//	| Subgraph) ";"?;
		@Override public ParserRule getRule() { return rule; }
		
		//(NodeStatement | EdgeStatement | AttributeStatement | Attribute | Subgraph) ";"?
		public Group getGroup() { return cGroup; }
		
		//NodeStatement | EdgeStatement | AttributeStatement | Attribute | Subgraph
		public Alternatives getAlternatives_0() { return cAlternatives_0; }
		
		//NodeStatement
		public RuleCall getNodeStatementParserRuleCall_0_0() { return cNodeStatementParserRuleCall_0_0; }
		
		//EdgeStatement
		public RuleCall getEdgeStatementParserRuleCall_0_1() { return cEdgeStatementParserRuleCall_0_1; }
		
		//AttributeStatement
		public RuleCall getAttributeStatementParserRuleCall_0_2() { return cAttributeStatementParserRuleCall_0_2; }
		
		//Attribute
		public RuleCall getAttributeParserRuleCall_0_3() { return cAttributeParserRuleCall_0_3; }
		
		//Subgraph
		public RuleCall getSubgraphParserRuleCall_0_4() { return cSubgraphParserRuleCall_0_4; }
		
		//";"?
		public Keyword getSemicolonKeyword_1() { return cSemicolonKeyword_1; }
	}
	public class AttributeElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.alg.graphviz.dot.GraphvizDot.Attribute");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Assignment cNameAssignment_0 = (Assignment)cGroup.eContents().get(0);
		private final RuleCall cNameDotIDParserRuleCall_0_0 = (RuleCall)cNameAssignment_0.eContents().get(0);
		private final Keyword cEqualsSignKeyword_1 = (Keyword)cGroup.eContents().get(1);
		private final Assignment cValueAssignment_2 = (Assignment)cGroup.eContents().get(2);
		private final RuleCall cValueDotIDParserRuleCall_2_0 = (RuleCall)cValueAssignment_2.eContents().get(0);
		
		//Attribute:
		//	name=DotID "=" value=DotID;
		@Override public ParserRule getRule() { return rule; }
		
		//name=DotID "=" value=DotID
		public Group getGroup() { return cGroup; }
		
		//name=DotID
		public Assignment getNameAssignment_0() { return cNameAssignment_0; }
		
		//DotID
		public RuleCall getNameDotIDParserRuleCall_0_0() { return cNameDotIDParserRuleCall_0_0; }
		
		//"="
		public Keyword getEqualsSignKeyword_1() { return cEqualsSignKeyword_1; }
		
		//value=DotID
		public Assignment getValueAssignment_2() { return cValueAssignment_2; }
		
		//DotID
		public RuleCall getValueDotIDParserRuleCall_2_0() { return cValueDotIDParserRuleCall_2_0; }
	}
	public class NodeStatementElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.alg.graphviz.dot.GraphvizDot.NodeStatement");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Assignment cNodeAssignment_0 = (Assignment)cGroup.eContents().get(0);
		private final RuleCall cNodeNodeParserRuleCall_0_0 = (RuleCall)cNodeAssignment_0.eContents().get(0);
		private final Group cGroup_1 = (Group)cGroup.eContents().get(1);
		private final Keyword cLeftSquareBracketKeyword_1_0 = (Keyword)cGroup_1.eContents().get(0);
		private final Group cGroup_1_1 = (Group)cGroup_1.eContents().get(1);
		private final Assignment cAttributesAssignment_1_1_0 = (Assignment)cGroup_1_1.eContents().get(0);
		private final RuleCall cAttributesListAttributeParserRuleCall_1_1_0_0 = (RuleCall)cAttributesAssignment_1_1_0.eContents().get(0);
		private final Group cGroup_1_1_1 = (Group)cGroup_1_1.eContents().get(1);
		private final Keyword cCommaKeyword_1_1_1_0 = (Keyword)cGroup_1_1_1.eContents().get(0);
		private final Assignment cAttributesAssignment_1_1_1_1 = (Assignment)cGroup_1_1_1.eContents().get(1);
		private final RuleCall cAttributesListAttributeParserRuleCall_1_1_1_1_0 = (RuleCall)cAttributesAssignment_1_1_1_1.eContents().get(0);
		private final Keyword cRightSquareBracketKeyword_1_2 = (Keyword)cGroup_1.eContents().get(2);
		
		///*~~~~Original definition~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// *  node_stmt:  node_id [ attr_list ]
		// */ NodeStatement:
		//	node=Node ("[" (attributes+=ListAttribute (","? attributes+=ListAttribute)*)? "]")?;
		@Override public ParserRule getRule() { return rule; }
		
		//node=Node ("[" (attributes+=ListAttribute (","? attributes+=ListAttribute)*)? "]")?
		public Group getGroup() { return cGroup; }
		
		//node=Node
		public Assignment getNodeAssignment_0() { return cNodeAssignment_0; }
		
		//Node
		public RuleCall getNodeNodeParserRuleCall_0_0() { return cNodeNodeParserRuleCall_0_0; }
		
		//("[" (attributes+=ListAttribute (","? attributes+=ListAttribute)*)? "]")?
		public Group getGroup_1() { return cGroup_1; }
		
		//"["
		public Keyword getLeftSquareBracketKeyword_1_0() { return cLeftSquareBracketKeyword_1_0; }
		
		//(attributes+=ListAttribute (","? attributes+=ListAttribute)*)?
		public Group getGroup_1_1() { return cGroup_1_1; }
		
		//attributes+=ListAttribute
		public Assignment getAttributesAssignment_1_1_0() { return cAttributesAssignment_1_1_0; }
		
		//ListAttribute
		public RuleCall getAttributesListAttributeParserRuleCall_1_1_0_0() { return cAttributesListAttributeParserRuleCall_1_1_0_0; }
		
		//(","? attributes+=ListAttribute)*
		public Group getGroup_1_1_1() { return cGroup_1_1_1; }
		
		//","?
		public Keyword getCommaKeyword_1_1_1_0() { return cCommaKeyword_1_1_1_0; }
		
		//attributes+=ListAttribute
		public Assignment getAttributesAssignment_1_1_1_1() { return cAttributesAssignment_1_1_1_1; }
		
		//ListAttribute
		public RuleCall getAttributesListAttributeParserRuleCall_1_1_1_1_0() { return cAttributesListAttributeParserRuleCall_1_1_1_1_0; }
		
		//"]"
		public Keyword getRightSquareBracketKeyword_1_2() { return cRightSquareBracketKeyword_1_2; }
	}
	public class NodeElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.alg.graphviz.dot.GraphvizDot.Node");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Assignment cNameAssignment_0 = (Assignment)cGroup.eContents().get(0);
		private final RuleCall cNameDotIDParserRuleCall_0_0 = (RuleCall)cNameAssignment_0.eContents().get(0);
		private final Group cGroup_1 = (Group)cGroup.eContents().get(1);
		private final Keyword cColonKeyword_1_0 = (Keyword)cGroup_1.eContents().get(0);
		private final Assignment cPortAssignment_1_1 = (Assignment)cGroup_1.eContents().get(1);
		private final RuleCall cPortPortParserRuleCall_1_1_0 = (RuleCall)cPortAssignment_1_1.eContents().get(0);
		
		///*~~~~Original definition~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// *  node_id:    ID [ port ]
		// */ Node:
		//	name=DotID (":" port=Port)?;
		@Override public ParserRule getRule() { return rule; }
		
		//name=DotID (":" port=Port)?
		public Group getGroup() { return cGroup; }
		
		//name=DotID
		public Assignment getNameAssignment_0() { return cNameAssignment_0; }
		
		//DotID
		public RuleCall getNameDotIDParserRuleCall_0_0() { return cNameDotIDParserRuleCall_0_0; }
		
		//(":" port=Port)?
		public Group getGroup_1() { return cGroup_1; }
		
		//":"
		public Keyword getColonKeyword_1_0() { return cColonKeyword_1_0; }
		
		//port=Port
		public Assignment getPortAssignment_1_1() { return cPortAssignment_1_1; }
		
		//Port
		public RuleCall getPortPortParserRuleCall_1_1_0() { return cPortPortParserRuleCall_1_1_0; }
	}
	public class EdgeStatementElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.alg.graphviz.dot.GraphvizDot.EdgeStatement");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Assignment cSourceNodeAssignment_0 = (Assignment)cGroup.eContents().get(0);
		private final RuleCall cSourceNodeNodeParserRuleCall_0_0 = (RuleCall)cSourceNodeAssignment_0.eContents().get(0);
		private final Assignment cEdgeTargetsAssignment_1 = (Assignment)cGroup.eContents().get(1);
		private final RuleCall cEdgeTargetsEdgeTargetParserRuleCall_1_0 = (RuleCall)cEdgeTargetsAssignment_1.eContents().get(0);
		private final Group cGroup_2 = (Group)cGroup.eContents().get(2);
		private final Keyword cLeftSquareBracketKeyword_2_0 = (Keyword)cGroup_2.eContents().get(0);
		private final Group cGroup_2_1 = (Group)cGroup_2.eContents().get(1);
		private final Assignment cAttributesAssignment_2_1_0 = (Assignment)cGroup_2_1.eContents().get(0);
		private final RuleCall cAttributesListAttributeParserRuleCall_2_1_0_0 = (RuleCall)cAttributesAssignment_2_1_0.eContents().get(0);
		private final Group cGroup_2_1_1 = (Group)cGroup_2_1.eContents().get(1);
		private final Keyword cCommaKeyword_2_1_1_0 = (Keyword)cGroup_2_1_1.eContents().get(0);
		private final Assignment cAttributesAssignment_2_1_1_1 = (Assignment)cGroup_2_1_1.eContents().get(1);
		private final RuleCall cAttributesListAttributeParserRuleCall_2_1_1_1_0 = (RuleCall)cAttributesAssignment_2_1_1_1.eContents().get(0);
		private final Keyword cRightSquareBracketKeyword_2_2 = (Keyword)cGroup_2.eContents().get(2);
		
		///*~~~~Original definition~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// *  edge_stmt:  (node_id | subgraph) edgeRHS [ attr_list ]
		// */ EdgeStatement:
		//	sourceNode=Node edgeTargets+=EdgeTarget+ ("[" (attributes+=ListAttribute (","? attributes+=ListAttribute)*)? "]")?;
		@Override public ParserRule getRule() { return rule; }
		
		//sourceNode=Node edgeTargets+=EdgeTarget+ ("[" (attributes+=ListAttribute (","? attributes+=ListAttribute)*)? "]")?
		public Group getGroup() { return cGroup; }
		
		//sourceNode=Node
		public Assignment getSourceNodeAssignment_0() { return cSourceNodeAssignment_0; }
		
		//Node
		public RuleCall getSourceNodeNodeParserRuleCall_0_0() { return cSourceNodeNodeParserRuleCall_0_0; }
		
		//edgeTargets+=EdgeTarget+
		public Assignment getEdgeTargetsAssignment_1() { return cEdgeTargetsAssignment_1; }
		
		//EdgeTarget
		public RuleCall getEdgeTargetsEdgeTargetParserRuleCall_1_0() { return cEdgeTargetsEdgeTargetParserRuleCall_1_0; }
		
		//("[" (attributes+=ListAttribute (","? attributes+=ListAttribute)*)? "]")?
		public Group getGroup_2() { return cGroup_2; }
		
		//"["
		public Keyword getLeftSquareBracketKeyword_2_0() { return cLeftSquareBracketKeyword_2_0; }
		
		//(attributes+=ListAttribute (","? attributes+=ListAttribute)*)?
		public Group getGroup_2_1() { return cGroup_2_1; }
		
		//attributes+=ListAttribute
		public Assignment getAttributesAssignment_2_1_0() { return cAttributesAssignment_2_1_0; }
		
		//ListAttribute
		public RuleCall getAttributesListAttributeParserRuleCall_2_1_0_0() { return cAttributesListAttributeParserRuleCall_2_1_0_0; }
		
		//(","? attributes+=ListAttribute)*
		public Group getGroup_2_1_1() { return cGroup_2_1_1; }
		
		//","?
		public Keyword getCommaKeyword_2_1_1_0() { return cCommaKeyword_2_1_1_0; }
		
		//attributes+=ListAttribute
		public Assignment getAttributesAssignment_2_1_1_1() { return cAttributesAssignment_2_1_1_1; }
		
		//ListAttribute
		public RuleCall getAttributesListAttributeParserRuleCall_2_1_1_1_0() { return cAttributesListAttributeParserRuleCall_2_1_1_1_0; }
		
		//"]"
		public Keyword getRightSquareBracketKeyword_2_2() { return cRightSquareBracketKeyword_2_2; }
	}
	public class EdgeTargetElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.alg.graphviz.dot.GraphvizDot.EdgeTarget");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Assignment cOperatorAssignment_0 = (Assignment)cGroup.eContents().get(0);
		private final RuleCall cOperatorEdgeOperatorEnumRuleCall_0_0 = (RuleCall)cOperatorAssignment_0.eContents().get(0);
		private final Alternatives cAlternatives_1 = (Alternatives)cGroup.eContents().get(1);
		private final Assignment cTargetSubgraphAssignment_1_0 = (Assignment)cAlternatives_1.eContents().get(0);
		private final RuleCall cTargetSubgraphSubgraphParserRuleCall_1_0_0 = (RuleCall)cTargetSubgraphAssignment_1_0.eContents().get(0);
		private final Assignment cTargetnodeAssignment_1_1 = (Assignment)cAlternatives_1.eContents().get(1);
		private final RuleCall cTargetnodeNodeParserRuleCall_1_1_0 = (RuleCall)cTargetnodeAssignment_1_1.eContents().get(0);
		
		///*~~~~Original definition~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// *  edgeRHS:    edgeop (node_id | subgraph) [ edgeRHS ]
		// */ EdgeTarget:
		//	operator=EdgeOperator (targetSubgraph=Subgraph | targetnode=Node);
		@Override public ParserRule getRule() { return rule; }
		
		//operator=EdgeOperator (targetSubgraph=Subgraph | targetnode=Node)
		public Group getGroup() { return cGroup; }
		
		//operator=EdgeOperator
		public Assignment getOperatorAssignment_0() { return cOperatorAssignment_0; }
		
		//EdgeOperator
		public RuleCall getOperatorEdgeOperatorEnumRuleCall_0_0() { return cOperatorEdgeOperatorEnumRuleCall_0_0; }
		
		//targetSubgraph=Subgraph | targetnode=Node
		public Alternatives getAlternatives_1() { return cAlternatives_1; }
		
		//targetSubgraph=Subgraph
		public Assignment getTargetSubgraphAssignment_1_0() { return cTargetSubgraphAssignment_1_0; }
		
		//Subgraph
		public RuleCall getTargetSubgraphSubgraphParserRuleCall_1_0_0() { return cTargetSubgraphSubgraphParserRuleCall_1_0_0; }
		
		//targetnode=Node
		public Assignment getTargetnodeAssignment_1_1() { return cTargetnodeAssignment_1_1; }
		
		//Node
		public RuleCall getTargetnodeNodeParserRuleCall_1_1_0() { return cTargetnodeNodeParserRuleCall_1_1_0; }
	}
	public class AttributeStatementElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.alg.graphviz.dot.GraphvizDot.AttributeStatement");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Assignment cTypeAssignment_0 = (Assignment)cGroup.eContents().get(0);
		private final RuleCall cTypeAttributeTypeEnumRuleCall_0_0 = (RuleCall)cTypeAssignment_0.eContents().get(0);
		private final Keyword cLeftSquareBracketKeyword_1 = (Keyword)cGroup.eContents().get(1);
		private final Group cGroup_2 = (Group)cGroup.eContents().get(2);
		private final Assignment cAttributesAssignment_2_0 = (Assignment)cGroup_2.eContents().get(0);
		private final RuleCall cAttributesListAttributeParserRuleCall_2_0_0 = (RuleCall)cAttributesAssignment_2_0.eContents().get(0);
		private final Group cGroup_2_1 = (Group)cGroup_2.eContents().get(1);
		private final Keyword cCommaKeyword_2_1_0 = (Keyword)cGroup_2_1.eContents().get(0);
		private final Assignment cAttributesAssignment_2_1_1 = (Assignment)cGroup_2_1.eContents().get(1);
		private final RuleCall cAttributesListAttributeParserRuleCall_2_1_1_0 = (RuleCall)cAttributesAssignment_2_1_1.eContents().get(0);
		private final Keyword cRightSquareBracketKeyword_3 = (Keyword)cGroup.eContents().get(3);
		
		///*~~~~Original definition~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// *  attr_stmt:  ('graph' | 'node' | 'edge') attr_list
		// */ AttributeStatement:
		//	type=AttributeType
		//	"[" (attributes+=ListAttribute (","? attributes+=ListAttribute)*)? "]";
		@Override public ParserRule getRule() { return rule; }
		
		//type=AttributeType "[" (attributes+=ListAttribute (","? attributes+=ListAttribute)*)? "]"
		public Group getGroup() { return cGroup; }
		
		//type=AttributeType
		public Assignment getTypeAssignment_0() { return cTypeAssignment_0; }
		
		//AttributeType
		public RuleCall getTypeAttributeTypeEnumRuleCall_0_0() { return cTypeAttributeTypeEnumRuleCall_0_0; }
		
		//"["
		public Keyword getLeftSquareBracketKeyword_1() { return cLeftSquareBracketKeyword_1; }
		
		//(attributes+=ListAttribute (","? attributes+=ListAttribute)*)?
		public Group getGroup_2() { return cGroup_2; }
		
		//attributes+=ListAttribute
		public Assignment getAttributesAssignment_2_0() { return cAttributesAssignment_2_0; }
		
		//ListAttribute
		public RuleCall getAttributesListAttributeParserRuleCall_2_0_0() { return cAttributesListAttributeParserRuleCall_2_0_0; }
		
		//(","? attributes+=ListAttribute)*
		public Group getGroup_2_1() { return cGroup_2_1; }
		
		//","?
		public Keyword getCommaKeyword_2_1_0() { return cCommaKeyword_2_1_0; }
		
		//attributes+=ListAttribute
		public Assignment getAttributesAssignment_2_1_1() { return cAttributesAssignment_2_1_1; }
		
		//ListAttribute
		public RuleCall getAttributesListAttributeParserRuleCall_2_1_1_0() { return cAttributesListAttributeParserRuleCall_2_1_1_0; }
		
		//"]"
		public Keyword getRightSquareBracketKeyword_3() { return cRightSquareBracketKeyword_3; }
	}
	public class SubgraphElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.alg.graphviz.dot.GraphvizDot.Subgraph");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Action cSubgraphAction_0 = (Action)cGroup.eContents().get(0);
		private final Group cGroup_1 = (Group)cGroup.eContents().get(1);
		private final Keyword cSubgraphKeyword_1_0 = (Keyword)cGroup_1.eContents().get(0);
		private final Assignment cNameAssignment_1_1 = (Assignment)cGroup_1.eContents().get(1);
		private final RuleCall cNameIDTerminalRuleCall_1_1_0 = (RuleCall)cNameAssignment_1_1.eContents().get(0);
		private final Keyword cLeftCurlyBracketKeyword_2 = (Keyword)cGroup.eContents().get(2);
		private final Assignment cStatementsAssignment_3 = (Assignment)cGroup.eContents().get(3);
		private final RuleCall cStatementsStatementParserRuleCall_3_0 = (RuleCall)cStatementsAssignment_3.eContents().get(0);
		private final Keyword cRightCurlyBracketKeyword_4 = (Keyword)cGroup.eContents().get(4);
		
		///*~~~~Original definition~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// *  subgraph:   [ 'subgraph' [ ID ] ] '{' stmt_list '}'
		// */ Subgraph:
		//	{Subgraph} ("subgraph" name=ID?)? "{" statements+=Statement* "}";
		@Override public ParserRule getRule() { return rule; }
		
		//{Subgraph} ("subgraph" name=ID?)? "{" statements+=Statement* "}"
		public Group getGroup() { return cGroup; }
		
		//{Subgraph}
		public Action getSubgraphAction_0() { return cSubgraphAction_0; }
		
		//("subgraph" name=ID?)?
		public Group getGroup_1() { return cGroup_1; }
		
		//"subgraph"
		public Keyword getSubgraphKeyword_1_0() { return cSubgraphKeyword_1_0; }
		
		//name=ID?
		public Assignment getNameAssignment_1_1() { return cNameAssignment_1_1; }
		
		//ID
		public RuleCall getNameIDTerminalRuleCall_1_1_0() { return cNameIDTerminalRuleCall_1_1_0; }
		
		//"{"
		public Keyword getLeftCurlyBracketKeyword_2() { return cLeftCurlyBracketKeyword_2; }
		
		//statements+=Statement*
		public Assignment getStatementsAssignment_3() { return cStatementsAssignment_3; }
		
		//Statement
		public RuleCall getStatementsStatementParserRuleCall_3_0() { return cStatementsStatementParserRuleCall_3_0; }
		
		//"}"
		public Keyword getRightCurlyBracketKeyword_4() { return cRightCurlyBracketKeyword_4; }
	}
	public class ListAttributeElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.alg.graphviz.dot.GraphvizDot.ListAttribute");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Assignment cNameAssignment_0 = (Assignment)cGroup.eContents().get(0);
		private final RuleCall cNameDotIDParserRuleCall_0_0 = (RuleCall)cNameAssignment_0.eContents().get(0);
		private final Group cGroup_1 = (Group)cGroup.eContents().get(1);
		private final Keyword cEqualsSignKeyword_1_0 = (Keyword)cGroup_1.eContents().get(0);
		private final Assignment cValueAssignment_1_1 = (Assignment)cGroup_1.eContents().get(1);
		private final RuleCall cValueDotIDParserRuleCall_1_1_0 = (RuleCall)cValueAssignment_1_1.eContents().get(0);
		
		///*~~~~Original definition~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// *  attr_list:  '[' [ a_list ] ']' [ attr_list ]
		// *  a_list:     ID [ '=' ID ] [ ',' ] [ a_list ]
		// */ ListAttribute Attribute:
		//	name=DotID ("=" value=DotID)?;
		@Override public ParserRule getRule() { return rule; }
		
		//name=DotID ("=" value=DotID)?
		public Group getGroup() { return cGroup; }
		
		//name=DotID
		public Assignment getNameAssignment_0() { return cNameAssignment_0; }
		
		//DotID
		public RuleCall getNameDotIDParserRuleCall_0_0() { return cNameDotIDParserRuleCall_0_0; }
		
		//("=" value=DotID)?
		public Group getGroup_1() { return cGroup_1; }
		
		//"="
		public Keyword getEqualsSignKeyword_1_0() { return cEqualsSignKeyword_1_0; }
		
		//value=DotID
		public Assignment getValueAssignment_1_1() { return cValueAssignment_1_1; }
		
		//DotID
		public RuleCall getValueDotIDParserRuleCall_1_1_0() { return cValueDotIDParserRuleCall_1_1_0; }
	}
	public class PortElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.alg.graphviz.dot.GraphvizDot.Port");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Assignment cNameAssignment_0 = (Assignment)cGroup.eContents().get(0);
		private final RuleCall cNameDotIDParserRuleCall_0_0 = (RuleCall)cNameAssignment_0.eContents().get(0);
		private final Group cGroup_1 = (Group)cGroup.eContents().get(1);
		private final Keyword cColonKeyword_1_0 = (Keyword)cGroup_1.eContents().get(0);
		private final Assignment cCompass_ptAssignment_1_1 = (Assignment)cGroup_1.eContents().get(1);
		private final RuleCall cCompass_ptIDTerminalRuleCall_1_1_0 = (RuleCall)cCompass_ptAssignment_1_1.eContents().get(0);
		
		///*~~~~Original definition~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// *  port:       ':' ID [ ':' compass_pt ] | ':' compass_pt
		// *  compass_pt: (n | ne | e | se | s | sw | w | nw | c | _)
		// */ Port:
		//	name=DotID (":" compass_pt=ID)?;
		@Override public ParserRule getRule() { return rule; }
		
		//name=DotID (":" compass_pt=ID)?
		public Group getGroup() { return cGroup; }
		
		//name=DotID
		public Assignment getNameAssignment_0() { return cNameAssignment_0; }
		
		//DotID
		public RuleCall getNameDotIDParserRuleCall_0_0() { return cNameDotIDParserRuleCall_0_0; }
		
		//(":" compass_pt=ID)?
		public Group getGroup_1() { return cGroup_1; }
		
		//":"
		public Keyword getColonKeyword_1_0() { return cColonKeyword_1_0; }
		
		//compass_pt=ID
		public Assignment getCompass_ptAssignment_1_1() { return cCompass_ptAssignment_1_1; }
		
		//ID
		public RuleCall getCompass_ptIDTerminalRuleCall_1_1_0() { return cCompass_ptIDTerminalRuleCall_1_1_0; }
	}
	public class DotIDElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.alg.graphviz.dot.GraphvizDot.DotID");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final RuleCall cIDTerminalRuleCall_0 = (RuleCall)cAlternatives.eContents().get(0);
		private final RuleCall cINTTerminalRuleCall_1 = (RuleCall)cAlternatives.eContents().get(1);
		private final RuleCall cFLOATTerminalRuleCall_2 = (RuleCall)cAlternatives.eContents().get(2);
		private final RuleCall cSTRINGTerminalRuleCall_3 = (RuleCall)cAlternatives.eContents().get(3);
		
		///*~~~~Original definition~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// *  An ID is one of the following:
		// *    - Any string of alphabetic ([a-zA-Z'200-'377]) characters, underscores
		// *      ('_') or digits ([0-9]), not beginning with a digit;
		// *    - a numeral [-]?(.[0-9]+ | [0-9]+(.[0-9]*)? );
		// *    - any double-quoted string ("...") possibly containing escaped quotes ('");
		// *    - an HTML string (<...>). 
		// */ DotID:
		//	ID | INT | FLOAT | STRING;
		@Override public ParserRule getRule() { return rule; }
		
		//ID | INT | FLOAT | STRING
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//ID
		public RuleCall getIDTerminalRuleCall_0() { return cIDTerminalRuleCall_0; }
		
		//INT
		public RuleCall getINTTerminalRuleCall_1() { return cINTTerminalRuleCall_1; }
		
		//FLOAT
		public RuleCall getFLOATTerminalRuleCall_2() { return cFLOATTerminalRuleCall_2; }
		
		//STRING
		public RuleCall getSTRINGTerminalRuleCall_3() { return cSTRINGTerminalRuleCall_3; }
	}
	
	public class EdgeOperatorElements extends AbstractEnumRuleElementFinder {
		private final EnumRule rule = (EnumRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.alg.graphviz.dot.GraphvizDot.EdgeOperator");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final EnumLiteralDeclaration cDirectedEnumLiteralDeclaration_0 = (EnumLiteralDeclaration)cAlternatives.eContents().get(0);
		private final Keyword cDirectedHyphenMinusGreaterThanSignKeyword_0_0 = (Keyword)cDirectedEnumLiteralDeclaration_0.eContents().get(0);
		private final EnumLiteralDeclaration cUndirectedEnumLiteralDeclaration_1 = (EnumLiteralDeclaration)cAlternatives.eContents().get(1);
		private final Keyword cUndirectedHyphenMinusHyphenMinusKeyword_1_0 = (Keyword)cUndirectedEnumLiteralDeclaration_1.eContents().get(0);
		
		//enum EdgeOperator:
		//	directed="->" | undirected="--";
		public EnumRule getRule() { return rule; }
		
		//directed="->" | undirected="--"
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//directed="->"
		public EnumLiteralDeclaration getDirectedEnumLiteralDeclaration_0() { return cDirectedEnumLiteralDeclaration_0; }
		
		//"->"
		public Keyword getDirectedHyphenMinusGreaterThanSignKeyword_0_0() { return cDirectedHyphenMinusGreaterThanSignKeyword_0_0; }
		
		//undirected="--"
		public EnumLiteralDeclaration getUndirectedEnumLiteralDeclaration_1() { return cUndirectedEnumLiteralDeclaration_1; }
		
		//"--"
		public Keyword getUndirectedHyphenMinusHyphenMinusKeyword_1_0() { return cUndirectedHyphenMinusHyphenMinusKeyword_1_0; }
	}
	public class GraphTypeElements extends AbstractEnumRuleElementFinder {
		private final EnumRule rule = (EnumRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.alg.graphviz.dot.GraphvizDot.GraphType");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final EnumLiteralDeclaration cGraphEnumLiteralDeclaration_0 = (EnumLiteralDeclaration)cAlternatives.eContents().get(0);
		private final Keyword cGraphGraphKeyword_0_0 = (Keyword)cGraphEnumLiteralDeclaration_0.eContents().get(0);
		private final EnumLiteralDeclaration cDigraphEnumLiteralDeclaration_1 = (EnumLiteralDeclaration)cAlternatives.eContents().get(1);
		private final Keyword cDigraphDigraphKeyword_1_0 = (Keyword)cDigraphEnumLiteralDeclaration_1.eContents().get(0);
		
		//enum GraphType:
		//	graph | digraph;
		public EnumRule getRule() { return rule; }
		
		//graph | digraph
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//graph
		public EnumLiteralDeclaration getGraphEnumLiteralDeclaration_0() { return cGraphEnumLiteralDeclaration_0; }
		
		//"graph"
		public Keyword getGraphGraphKeyword_0_0() { return cGraphGraphKeyword_0_0; }
		
		//digraph
		public EnumLiteralDeclaration getDigraphEnumLiteralDeclaration_1() { return cDigraphEnumLiteralDeclaration_1; }
		
		//"digraph"
		public Keyword getDigraphDigraphKeyword_1_0() { return cDigraphDigraphKeyword_1_0; }
	}
	public class AttributeTypeElements extends AbstractEnumRuleElementFinder {
		private final EnumRule rule = (EnumRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.alg.graphviz.dot.GraphvizDot.AttributeType");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final EnumLiteralDeclaration cGraphEnumLiteralDeclaration_0 = (EnumLiteralDeclaration)cAlternatives.eContents().get(0);
		private final Keyword cGraphGraphKeyword_0_0 = (Keyword)cGraphEnumLiteralDeclaration_0.eContents().get(0);
		private final EnumLiteralDeclaration cNodeEnumLiteralDeclaration_1 = (EnumLiteralDeclaration)cAlternatives.eContents().get(1);
		private final Keyword cNodeNodeKeyword_1_0 = (Keyword)cNodeEnumLiteralDeclaration_1.eContents().get(0);
		private final EnumLiteralDeclaration cEdgeEnumLiteralDeclaration_2 = (EnumLiteralDeclaration)cAlternatives.eContents().get(2);
		private final Keyword cEdgeEdgeKeyword_2_0 = (Keyword)cEdgeEnumLiteralDeclaration_2.eContents().get(0);
		
		//enum AttributeType:
		//	graph | node | edge;
		public EnumRule getRule() { return rule; }
		
		//graph | node | edge
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//graph
		public EnumLiteralDeclaration getGraphEnumLiteralDeclaration_0() { return cGraphEnumLiteralDeclaration_0; }
		
		//"graph"
		public Keyword getGraphGraphKeyword_0_0() { return cGraphGraphKeyword_0_0; }
		
		//node
		public EnumLiteralDeclaration getNodeEnumLiteralDeclaration_1() { return cNodeEnumLiteralDeclaration_1; }
		
		//"node"
		public Keyword getNodeNodeKeyword_1_0() { return cNodeNodeKeyword_1_0; }
		
		//edge
		public EnumLiteralDeclaration getEdgeEnumLiteralDeclaration_2() { return cEdgeEnumLiteralDeclaration_2; }
		
		//"edge"
		public Keyword getEdgeEdgeKeyword_2_0() { return cEdgeEdgeKeyword_2_0; }
	}
	
	private final GraphvizModelElements pGraphvizModel;
	private final GraphElements pGraph;
	private final StatementElements pStatement;
	private final AttributeElements pAttribute;
	private final NodeStatementElements pNodeStatement;
	private final NodeElements pNode;
	private final EdgeStatementElements pEdgeStatement;
	private final EdgeTargetElements pEdgeTarget;
	private final AttributeStatementElements pAttributeStatement;
	private final SubgraphElements pSubgraph;
	private final ListAttributeElements pListAttribute;
	private final PortElements pPort;
	private final DotIDElements pDotID;
	private final EdgeOperatorElements eEdgeOperator;
	private final GraphTypeElements eGraphType;
	private final AttributeTypeElements eAttributeType;
	private final TerminalRule tINT;
	private final TerminalRule tID;
	private final TerminalRule tFLOAT;
	private final TerminalRule tSTRING;
	private final TerminalRule tPREC_LINE;
	
	private final Grammar grammar;
	
	private final TerminalsGrammarAccess gaTerminals;

	@Inject
	public GraphvizDotGrammarAccess(GrammarProvider grammarProvider,
			TerminalsGrammarAccess gaTerminals) {
		this.grammar = internalFindGrammar(grammarProvider);
		this.gaTerminals = gaTerminals;
		this.pGraphvizModel = new GraphvizModelElements();
		this.pGraph = new GraphElements();
		this.pStatement = new StatementElements();
		this.pAttribute = new AttributeElements();
		this.pNodeStatement = new NodeStatementElements();
		this.pNode = new NodeElements();
		this.pEdgeStatement = new EdgeStatementElements();
		this.pEdgeTarget = new EdgeTargetElements();
		this.pAttributeStatement = new AttributeStatementElements();
		this.pSubgraph = new SubgraphElements();
		this.pListAttribute = new ListAttributeElements();
		this.pPort = new PortElements();
		this.pDotID = new DotIDElements();
		this.eEdgeOperator = new EdgeOperatorElements();
		this.eGraphType = new GraphTypeElements();
		this.eAttributeType = new AttributeTypeElements();
		this.tINT = (TerminalRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.alg.graphviz.dot.GraphvizDot.INT");
		this.tID = (TerminalRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.alg.graphviz.dot.GraphvizDot.ID");
		this.tFLOAT = (TerminalRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.alg.graphviz.dot.GraphvizDot.FLOAT");
		this.tSTRING = (TerminalRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.alg.graphviz.dot.GraphvizDot.STRING");
		this.tPREC_LINE = (TerminalRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.alg.graphviz.dot.GraphvizDot.PREC_LINE");
	}
	
	protected Grammar internalFindGrammar(GrammarProvider grammarProvider) {
		Grammar grammar = grammarProvider.getGrammar(this);
		while (grammar != null) {
			if ("org.eclipse.elk.alg.graphviz.dot.GraphvizDot".equals(grammar.getName())) {
				return grammar;
			}
			List<Grammar> grammars = grammar.getUsedGrammars();
			if (!grammars.isEmpty()) {
				grammar = grammars.iterator().next();
			} else {
				return null;
			}
		}
		return grammar;
	}
	
	@Override
	public Grammar getGrammar() {
		return grammar;
	}
	
	
	public TerminalsGrammarAccess getTerminalsGrammarAccess() {
		return gaTerminals;
	}

	
	//GraphvizModel:
	//	graphs+=Graph*;
	public GraphvizModelElements getGraphvizModelAccess() {
		return pGraphvizModel;
	}
	
	public ParserRule getGraphvizModelRule() {
		return getGraphvizModelAccess().getRule();
	}
	
	///*~~~~Original definition~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// *  graph:      [ 'strict' ] ('graph' | 'digraph') [ ID ] '{' stmt_list '}'
	// */ Graph:
	//	strict?="strict"? type=GraphType name=DotID?
	//	"{" statements+=Statement* "}";
	public GraphElements getGraphAccess() {
		return pGraph;
	}
	
	public ParserRule getGraphRule() {
		return getGraphAccess().getRule();
	}
	
	///*~~~~Original definition~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// *  stmt_list:  [ stmt [ ';' ] [ stmt_list ] ]
	// *  stmt:       node_stmt | edge_stmt | attr_stmt | ID '=' ID | subgraph
	// */ Statement:
	//	(NodeStatement | EdgeStatement | AttributeStatement | Attribute
	//	| Subgraph) ";"?;
	public StatementElements getStatementAccess() {
		return pStatement;
	}
	
	public ParserRule getStatementRule() {
		return getStatementAccess().getRule();
	}
	
	//Attribute:
	//	name=DotID "=" value=DotID;
	public AttributeElements getAttributeAccess() {
		return pAttribute;
	}
	
	public ParserRule getAttributeRule() {
		return getAttributeAccess().getRule();
	}
	
	///*~~~~Original definition~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// *  node_stmt:  node_id [ attr_list ]
	// */ NodeStatement:
	//	node=Node ("[" (attributes+=ListAttribute (","? attributes+=ListAttribute)*)? "]")?;
	public NodeStatementElements getNodeStatementAccess() {
		return pNodeStatement;
	}
	
	public ParserRule getNodeStatementRule() {
		return getNodeStatementAccess().getRule();
	}
	
	///*~~~~Original definition~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// *  node_id:    ID [ port ]
	// */ Node:
	//	name=DotID (":" port=Port)?;
	public NodeElements getNodeAccess() {
		return pNode;
	}
	
	public ParserRule getNodeRule() {
		return getNodeAccess().getRule();
	}
	
	///*~~~~Original definition~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// *  edge_stmt:  (node_id | subgraph) edgeRHS [ attr_list ]
	// */ EdgeStatement:
	//	sourceNode=Node edgeTargets+=EdgeTarget+ ("[" (attributes+=ListAttribute (","? attributes+=ListAttribute)*)? "]")?;
	public EdgeStatementElements getEdgeStatementAccess() {
		return pEdgeStatement;
	}
	
	public ParserRule getEdgeStatementRule() {
		return getEdgeStatementAccess().getRule();
	}
	
	///*~~~~Original definition~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// *  edgeRHS:    edgeop (node_id | subgraph) [ edgeRHS ]
	// */ EdgeTarget:
	//	operator=EdgeOperator (targetSubgraph=Subgraph | targetnode=Node);
	public EdgeTargetElements getEdgeTargetAccess() {
		return pEdgeTarget;
	}
	
	public ParserRule getEdgeTargetRule() {
		return getEdgeTargetAccess().getRule();
	}
	
	///*~~~~Original definition~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// *  attr_stmt:  ('graph' | 'node' | 'edge') attr_list
	// */ AttributeStatement:
	//	type=AttributeType
	//	"[" (attributes+=ListAttribute (","? attributes+=ListAttribute)*)? "]";
	public AttributeStatementElements getAttributeStatementAccess() {
		return pAttributeStatement;
	}
	
	public ParserRule getAttributeStatementRule() {
		return getAttributeStatementAccess().getRule();
	}
	
	///*~~~~Original definition~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// *  subgraph:   [ 'subgraph' [ ID ] ] '{' stmt_list '}'
	// */ Subgraph:
	//	{Subgraph} ("subgraph" name=ID?)? "{" statements+=Statement* "}";
	public SubgraphElements getSubgraphAccess() {
		return pSubgraph;
	}
	
	public ParserRule getSubgraphRule() {
		return getSubgraphAccess().getRule();
	}
	
	///*~~~~Original definition~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// *  attr_list:  '[' [ a_list ] ']' [ attr_list ]
	// *  a_list:     ID [ '=' ID ] [ ',' ] [ a_list ]
	// */ ListAttribute Attribute:
	//	name=DotID ("=" value=DotID)?;
	public ListAttributeElements getListAttributeAccess() {
		return pListAttribute;
	}
	
	public ParserRule getListAttributeRule() {
		return getListAttributeAccess().getRule();
	}
	
	///*~~~~Original definition~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// *  port:       ':' ID [ ':' compass_pt ] | ':' compass_pt
	// *  compass_pt: (n | ne | e | se | s | sw | w | nw | c | _)
	// */ Port:
	//	name=DotID (":" compass_pt=ID)?;
	public PortElements getPortAccess() {
		return pPort;
	}
	
	public ParserRule getPortRule() {
		return getPortAccess().getRule();
	}
	
	///*~~~~Original definition~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// *  An ID is one of the following:
	// *    - Any string of alphabetic ([a-zA-Z'200-'377]) characters, underscores
	// *      ('_') or digits ([0-9]), not beginning with a digit;
	// *    - a numeral [-]?(.[0-9]+ | [0-9]+(.[0-9]*)? );
	// *    - any double-quoted string ("...") possibly containing escaped quotes ('");
	// *    - an HTML string (<...>). 
	// */ DotID:
	//	ID | INT | FLOAT | STRING;
	public DotIDElements getDotIDAccess() {
		return pDotID;
	}
	
	public ParserRule getDotIDRule() {
		return getDotIDAccess().getRule();
	}
	
	//enum EdgeOperator:
	//	directed="->" | undirected="--";
	public EdgeOperatorElements getEdgeOperatorAccess() {
		return eEdgeOperator;
	}
	
	public EnumRule getEdgeOperatorRule() {
		return getEdgeOperatorAccess().getRule();
	}
	
	//enum GraphType:
	//	graph | digraph;
	public GraphTypeElements getGraphTypeAccess() {
		return eGraphType;
	}
	
	public EnumRule getGraphTypeRule() {
		return getGraphTypeAccess().getRule();
	}
	
	//enum AttributeType:
	//	graph | node | edge;
	public AttributeTypeElements getAttributeTypeAccess() {
		return eAttributeType;
	}
	
	public EnumRule getAttributeTypeRule() {
		return getAttributeTypeAccess().getRule();
	}
	
	//terminal INT returns ecore::EInt:
	//	'-'? '0'..'9'+;
	public TerminalRule getINTRule() {
		return tINT;
	}
	
	//terminal ID:
	//	('a'..'z' | 'A'..'Z' | '_') ('a'..'z' | 'A'..'Z' | '_' | '0'..'9')*;
	public TerminalRule getIDRule() {
		return tID;
	}
	
	//terminal FLOAT:
	//	'-'? ('0'..'9'* '.' '0'..'9'+);
	public TerminalRule getFLOATRule() {
		return tFLOAT;
	}
	
	//terminal STRING:
	//	'"' ('\\' '"' | !'"')* '"';
	public TerminalRule getSTRINGRule() {
		return tSTRING;
	}
	
	//terminal PREC_LINE:
	//	'#' !('\n' | '\r')* ('\r'? '\n')?;
	public TerminalRule getPREC_LINERule() {
		return tPREC_LINE;
	}
	
	//terminal ML_COMMENT:
	//	'/*'->'*/';
	public TerminalRule getML_COMMENTRule() {
		return gaTerminals.getML_COMMENTRule();
	}
	
	//terminal SL_COMMENT:
	//	'//' !('\n' | '\r')* ('\r'? '\n')?;
	public TerminalRule getSL_COMMENTRule() {
		return gaTerminals.getSL_COMMENTRule();
	}
	
	//terminal WS:
	//	' ' | '\t' | '\r' | '\n'+;
	public TerminalRule getWSRule() {
		return gaTerminals.getWSRule();
	}
	
	//terminal ANY_OTHER:
	//	.;
	public TerminalRule getANY_OTHERRule() {
		return gaTerminals.getANY_OTHERRule();
	}
}
