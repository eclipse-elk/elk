/*******************************************************************************
 * Copyright (c) 2016 TypeFox GmbH (http://www.typefox.io) and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.graph.text.services;

import com.google.inject.Inject;
import com.google.inject.Singleton;
import java.util.List;
import org.eclipse.xtext.Action;
import org.eclipse.xtext.Alternatives;
import org.eclipse.xtext.Assignment;
import org.eclipse.xtext.CrossReference;
import org.eclipse.xtext.Grammar;
import org.eclipse.xtext.GrammarUtil;
import org.eclipse.xtext.Group;
import org.eclipse.xtext.Keyword;
import org.eclipse.xtext.ParserRule;
import org.eclipse.xtext.RuleCall;
import org.eclipse.xtext.TerminalRule;
import org.eclipse.xtext.UnorderedGroup;
import org.eclipse.xtext.common.services.TerminalsGrammarAccess;
import org.eclipse.xtext.service.AbstractElementFinder.AbstractGrammarElementFinder;
import org.eclipse.xtext.service.GrammarProvider;

@Singleton
public class ElkGraphGrammarAccess extends AbstractGrammarElementFinder {
	
	public class RootNodeElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.text.ElkGraph.RootNode");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Action cElkNodeAction_0 = (Action)cGroup.eContents().get(0);
		private final Group cGroup_1 = (Group)cGroup.eContents().get(1);
		private final Keyword cGraphKeyword_1_0 = (Keyword)cGroup_1.eContents().get(0);
		private final Assignment cIdentifierAssignment_1_1 = (Assignment)cGroup_1.eContents().get(1);
		private final RuleCall cIdentifierIDTerminalRuleCall_1_1_0 = (RuleCall)cIdentifierAssignment_1_1.eContents().get(0);
		private final RuleCall cShapeLayoutParserRuleCall_2 = (RuleCall)cGroup.eContents().get(2);
		private final Assignment cPropertiesAssignment_3 = (Assignment)cGroup.eContents().get(3);
		private final RuleCall cPropertiesPropertyParserRuleCall_3_0 = (RuleCall)cPropertiesAssignment_3.eContents().get(0);
		private final Alternatives cAlternatives_4 = (Alternatives)cGroup.eContents().get(4);
		private final Assignment cLabelsAssignment_4_0 = (Assignment)cAlternatives_4.eContents().get(0);
		private final RuleCall cLabelsElkLabelParserRuleCall_4_0_0 = (RuleCall)cLabelsAssignment_4_0.eContents().get(0);
		private final Assignment cPortsAssignment_4_1 = (Assignment)cAlternatives_4.eContents().get(1);
		private final RuleCall cPortsElkPortParserRuleCall_4_1_0 = (RuleCall)cPortsAssignment_4_1.eContents().get(0);
		private final Assignment cChildrenAssignment_4_2 = (Assignment)cAlternatives_4.eContents().get(2);
		private final RuleCall cChildrenElkNodeParserRuleCall_4_2_0 = (RuleCall)cChildrenAssignment_4_2.eContents().get(0);
		private final Assignment cContainedEdgesAssignment_4_3 = (Assignment)cAlternatives_4.eContents().get(3);
		private final RuleCall cContainedEdgesElkEdgeParserRuleCall_4_3_0 = (RuleCall)cContainedEdgesAssignment_4_3.eContents().get(0);
		
		//RootNode ElkNode:
		//	{ElkNode} ('graph' identifier=ID)?
		//	ShapeLayout?
		//	properties+=Property* (labels+=ElkLabel | ports+=ElkPort | children+=ElkNode | containedEdges+=ElkEdge)*;
		@Override public ParserRule getRule() { return rule; }
		
		//{ElkNode} ('graph' identifier=ID)? ShapeLayout? properties+=Property* (labels+=ElkLabel | ports+=ElkPort |
		//children+=ElkNode | containedEdges+=ElkEdge)*
		public Group getGroup() { return cGroup; }
		
		//{ElkNode}
		public Action getElkNodeAction_0() { return cElkNodeAction_0; }
		
		//('graph' identifier=ID)?
		public Group getGroup_1() { return cGroup_1; }
		
		//'graph'
		public Keyword getGraphKeyword_1_0() { return cGraphKeyword_1_0; }
		
		//identifier=ID
		public Assignment getIdentifierAssignment_1_1() { return cIdentifierAssignment_1_1; }
		
		//ID
		public RuleCall getIdentifierIDTerminalRuleCall_1_1_0() { return cIdentifierIDTerminalRuleCall_1_1_0; }
		
		//ShapeLayout?
		public RuleCall getShapeLayoutParserRuleCall_2() { return cShapeLayoutParserRuleCall_2; }
		
		//properties+=Property*
		public Assignment getPropertiesAssignment_3() { return cPropertiesAssignment_3; }
		
		//Property
		public RuleCall getPropertiesPropertyParserRuleCall_3_0() { return cPropertiesPropertyParserRuleCall_3_0; }
		
		//(labels+=ElkLabel | ports+=ElkPort | children+=ElkNode | containedEdges+=ElkEdge)*
		public Alternatives getAlternatives_4() { return cAlternatives_4; }
		
		//labels+=ElkLabel
		public Assignment getLabelsAssignment_4_0() { return cLabelsAssignment_4_0; }
		
		//ElkLabel
		public RuleCall getLabelsElkLabelParserRuleCall_4_0_0() { return cLabelsElkLabelParserRuleCall_4_0_0; }
		
		//ports+=ElkPort
		public Assignment getPortsAssignment_4_1() { return cPortsAssignment_4_1; }
		
		//ElkPort
		public RuleCall getPortsElkPortParserRuleCall_4_1_0() { return cPortsElkPortParserRuleCall_4_1_0; }
		
		//children+=ElkNode
		public Assignment getChildrenAssignment_4_2() { return cChildrenAssignment_4_2; }
		
		//ElkNode
		public RuleCall getChildrenElkNodeParserRuleCall_4_2_0() { return cChildrenElkNodeParserRuleCall_4_2_0; }
		
		//containedEdges+=ElkEdge
		public Assignment getContainedEdgesAssignment_4_3() { return cContainedEdgesAssignment_4_3; }
		
		//ElkEdge
		public RuleCall getContainedEdgesElkEdgeParserRuleCall_4_3_0() { return cContainedEdgesElkEdgeParserRuleCall_4_3_0; }
	}
	public class ElkNodeElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.text.ElkGraph.ElkNode");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Keyword cNodeKeyword_0 = (Keyword)cGroup.eContents().get(0);
		private final Assignment cIdentifierAssignment_1 = (Assignment)cGroup.eContents().get(1);
		private final RuleCall cIdentifierIDTerminalRuleCall_1_0 = (RuleCall)cIdentifierAssignment_1.eContents().get(0);
		private final Group cGroup_2 = (Group)cGroup.eContents().get(2);
		private final Keyword cLeftCurlyBracketKeyword_2_0 = (Keyword)cGroup_2.eContents().get(0);
		private final RuleCall cShapeLayoutParserRuleCall_2_1 = (RuleCall)cGroup_2.eContents().get(1);
		private final Assignment cPropertiesAssignment_2_2 = (Assignment)cGroup_2.eContents().get(2);
		private final RuleCall cPropertiesPropertyParserRuleCall_2_2_0 = (RuleCall)cPropertiesAssignment_2_2.eContents().get(0);
		private final Alternatives cAlternatives_2_3 = (Alternatives)cGroup_2.eContents().get(3);
		private final Assignment cLabelsAssignment_2_3_0 = (Assignment)cAlternatives_2_3.eContents().get(0);
		private final RuleCall cLabelsElkLabelParserRuleCall_2_3_0_0 = (RuleCall)cLabelsAssignment_2_3_0.eContents().get(0);
		private final Assignment cPortsAssignment_2_3_1 = (Assignment)cAlternatives_2_3.eContents().get(1);
		private final RuleCall cPortsElkPortParserRuleCall_2_3_1_0 = (RuleCall)cPortsAssignment_2_3_1.eContents().get(0);
		private final Assignment cChildrenAssignment_2_3_2 = (Assignment)cAlternatives_2_3.eContents().get(2);
		private final RuleCall cChildrenElkNodeParserRuleCall_2_3_2_0 = (RuleCall)cChildrenAssignment_2_3_2.eContents().get(0);
		private final Assignment cContainedEdgesAssignment_2_3_3 = (Assignment)cAlternatives_2_3.eContents().get(3);
		private final RuleCall cContainedEdgesElkEdgeParserRuleCall_2_3_3_0 = (RuleCall)cContainedEdgesAssignment_2_3_3.eContents().get(0);
		private final Keyword cRightCurlyBracketKeyword_2_4 = (Keyword)cGroup_2.eContents().get(4);
		
		//ElkNode:
		//	'node' identifier=ID ('{'
		//	ShapeLayout?
		//	properties+=Property* (labels+=ElkLabel | ports+=ElkPort | children+=ElkNode | containedEdges+=ElkEdge)*
		//	'}')?;
		@Override public ParserRule getRule() { return rule; }
		
		//'node' identifier=ID ('{' ShapeLayout? properties+=Property* (labels+=ElkLabel | ports+=ElkPort | children+=ElkNode |
		//containedEdges+=ElkEdge)* '}')?
		public Group getGroup() { return cGroup; }
		
		//'node'
		public Keyword getNodeKeyword_0() { return cNodeKeyword_0; }
		
		//identifier=ID
		public Assignment getIdentifierAssignment_1() { return cIdentifierAssignment_1; }
		
		//ID
		public RuleCall getIdentifierIDTerminalRuleCall_1_0() { return cIdentifierIDTerminalRuleCall_1_0; }
		
		//('{' ShapeLayout? properties+=Property* (labels+=ElkLabel | ports+=ElkPort | children+=ElkNode |
		//containedEdges+=ElkEdge)* '}')?
		public Group getGroup_2() { return cGroup_2; }
		
		//'{'
		public Keyword getLeftCurlyBracketKeyword_2_0() { return cLeftCurlyBracketKeyword_2_0; }
		
		//ShapeLayout?
		public RuleCall getShapeLayoutParserRuleCall_2_1() { return cShapeLayoutParserRuleCall_2_1; }
		
		//properties+=Property*
		public Assignment getPropertiesAssignment_2_2() { return cPropertiesAssignment_2_2; }
		
		//Property
		public RuleCall getPropertiesPropertyParserRuleCall_2_2_0() { return cPropertiesPropertyParserRuleCall_2_2_0; }
		
		//(labels+=ElkLabel | ports+=ElkPort | children+=ElkNode | containedEdges+=ElkEdge)*
		public Alternatives getAlternatives_2_3() { return cAlternatives_2_3; }
		
		//labels+=ElkLabel
		public Assignment getLabelsAssignment_2_3_0() { return cLabelsAssignment_2_3_0; }
		
		//ElkLabel
		public RuleCall getLabelsElkLabelParserRuleCall_2_3_0_0() { return cLabelsElkLabelParserRuleCall_2_3_0_0; }
		
		//ports+=ElkPort
		public Assignment getPortsAssignment_2_3_1() { return cPortsAssignment_2_3_1; }
		
		//ElkPort
		public RuleCall getPortsElkPortParserRuleCall_2_3_1_0() { return cPortsElkPortParserRuleCall_2_3_1_0; }
		
		//children+=ElkNode
		public Assignment getChildrenAssignment_2_3_2() { return cChildrenAssignment_2_3_2; }
		
		//ElkNode
		public RuleCall getChildrenElkNodeParserRuleCall_2_3_2_0() { return cChildrenElkNodeParserRuleCall_2_3_2_0; }
		
		//containedEdges+=ElkEdge
		public Assignment getContainedEdgesAssignment_2_3_3() { return cContainedEdgesAssignment_2_3_3; }
		
		//ElkEdge
		public RuleCall getContainedEdgesElkEdgeParserRuleCall_2_3_3_0() { return cContainedEdgesElkEdgeParserRuleCall_2_3_3_0; }
		
		//'}'
		public Keyword getRightCurlyBracketKeyword_2_4() { return cRightCurlyBracketKeyword_2_4; }
	}
	public class ElkLabelElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.text.ElkGraph.ElkLabel");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Keyword cLabelKeyword_0 = (Keyword)cGroup.eContents().get(0);
		private final Group cGroup_1 = (Group)cGroup.eContents().get(1);
		private final Assignment cIdentifierAssignment_1_0 = (Assignment)cGroup_1.eContents().get(0);
		private final RuleCall cIdentifierIDTerminalRuleCall_1_0_0 = (RuleCall)cIdentifierAssignment_1_0.eContents().get(0);
		private final Keyword cColonKeyword_1_1 = (Keyword)cGroup_1.eContents().get(1);
		private final Assignment cTextAssignment_2 = (Assignment)cGroup.eContents().get(2);
		private final RuleCall cTextSTRINGTerminalRuleCall_2_0 = (RuleCall)cTextAssignment_2.eContents().get(0);
		private final Group cGroup_3 = (Group)cGroup.eContents().get(3);
		private final Keyword cLeftCurlyBracketKeyword_3_0 = (Keyword)cGroup_3.eContents().get(0);
		private final RuleCall cShapeLayoutParserRuleCall_3_1 = (RuleCall)cGroup_3.eContents().get(1);
		private final Assignment cPropertiesAssignment_3_2 = (Assignment)cGroup_3.eContents().get(2);
		private final RuleCall cPropertiesPropertyParserRuleCall_3_2_0 = (RuleCall)cPropertiesAssignment_3_2.eContents().get(0);
		private final Assignment cLabelsAssignment_3_3 = (Assignment)cGroup_3.eContents().get(3);
		private final RuleCall cLabelsElkLabelParserRuleCall_3_3_0 = (RuleCall)cLabelsAssignment_3_3.eContents().get(0);
		private final Keyword cRightCurlyBracketKeyword_3_4 = (Keyword)cGroup_3.eContents().get(4);
		
		//ElkLabel:
		//	'label' (identifier=ID ':')? text=STRING ('{'
		//	ShapeLayout?
		//	properties+=Property*
		//	labels+=ElkLabel*
		//	'}')?;
		@Override public ParserRule getRule() { return rule; }
		
		//'label' (identifier=ID ':')? text=STRING ('{' ShapeLayout? properties+=Property* labels+=ElkLabel* '}')?
		public Group getGroup() { return cGroup; }
		
		//'label'
		public Keyword getLabelKeyword_0() { return cLabelKeyword_0; }
		
		//(identifier=ID ':')?
		public Group getGroup_1() { return cGroup_1; }
		
		//identifier=ID
		public Assignment getIdentifierAssignment_1_0() { return cIdentifierAssignment_1_0; }
		
		//ID
		public RuleCall getIdentifierIDTerminalRuleCall_1_0_0() { return cIdentifierIDTerminalRuleCall_1_0_0; }
		
		//':'
		public Keyword getColonKeyword_1_1() { return cColonKeyword_1_1; }
		
		//text=STRING
		public Assignment getTextAssignment_2() { return cTextAssignment_2; }
		
		//STRING
		public RuleCall getTextSTRINGTerminalRuleCall_2_0() { return cTextSTRINGTerminalRuleCall_2_0; }
		
		//('{' ShapeLayout? properties+=Property* labels+=ElkLabel* '}')?
		public Group getGroup_3() { return cGroup_3; }
		
		//'{'
		public Keyword getLeftCurlyBracketKeyword_3_0() { return cLeftCurlyBracketKeyword_3_0; }
		
		//ShapeLayout?
		public RuleCall getShapeLayoutParserRuleCall_3_1() { return cShapeLayoutParserRuleCall_3_1; }
		
		//properties+=Property*
		public Assignment getPropertiesAssignment_3_2() { return cPropertiesAssignment_3_2; }
		
		//Property
		public RuleCall getPropertiesPropertyParserRuleCall_3_2_0() { return cPropertiesPropertyParserRuleCall_3_2_0; }
		
		//labels+=ElkLabel*
		public Assignment getLabelsAssignment_3_3() { return cLabelsAssignment_3_3; }
		
		//ElkLabel
		public RuleCall getLabelsElkLabelParserRuleCall_3_3_0() { return cLabelsElkLabelParserRuleCall_3_3_0; }
		
		//'}'
		public Keyword getRightCurlyBracketKeyword_3_4() { return cRightCurlyBracketKeyword_3_4; }
	}
	public class ElkPortElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.text.ElkGraph.ElkPort");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Keyword cPortKeyword_0 = (Keyword)cGroup.eContents().get(0);
		private final Assignment cIdentifierAssignment_1 = (Assignment)cGroup.eContents().get(1);
		private final RuleCall cIdentifierIDTerminalRuleCall_1_0 = (RuleCall)cIdentifierAssignment_1.eContents().get(0);
		private final Group cGroup_2 = (Group)cGroup.eContents().get(2);
		private final Keyword cLeftCurlyBracketKeyword_2_0 = (Keyword)cGroup_2.eContents().get(0);
		private final RuleCall cShapeLayoutParserRuleCall_2_1 = (RuleCall)cGroup_2.eContents().get(1);
		private final Assignment cPropertiesAssignment_2_2 = (Assignment)cGroup_2.eContents().get(2);
		private final RuleCall cPropertiesPropertyParserRuleCall_2_2_0 = (RuleCall)cPropertiesAssignment_2_2.eContents().get(0);
		private final Assignment cLabelsAssignment_2_3 = (Assignment)cGroup_2.eContents().get(3);
		private final RuleCall cLabelsElkLabelParserRuleCall_2_3_0 = (RuleCall)cLabelsAssignment_2_3.eContents().get(0);
		private final Keyword cRightCurlyBracketKeyword_2_4 = (Keyword)cGroup_2.eContents().get(4);
		
		//ElkPort:
		//	'port' identifier=ID ('{'
		//	ShapeLayout?
		//	properties+=Property*
		//	labels+=ElkLabel*
		//	'}')?;
		@Override public ParserRule getRule() { return rule; }
		
		//'port' identifier=ID ('{' ShapeLayout? properties+=Property* labels+=ElkLabel* '}')?
		public Group getGroup() { return cGroup; }
		
		//'port'
		public Keyword getPortKeyword_0() { return cPortKeyword_0; }
		
		//identifier=ID
		public Assignment getIdentifierAssignment_1() { return cIdentifierAssignment_1; }
		
		//ID
		public RuleCall getIdentifierIDTerminalRuleCall_1_0() { return cIdentifierIDTerminalRuleCall_1_0; }
		
		//('{' ShapeLayout? properties+=Property* labels+=ElkLabel* '}')?
		public Group getGroup_2() { return cGroup_2; }
		
		//'{'
		public Keyword getLeftCurlyBracketKeyword_2_0() { return cLeftCurlyBracketKeyword_2_0; }
		
		//ShapeLayout?
		public RuleCall getShapeLayoutParserRuleCall_2_1() { return cShapeLayoutParserRuleCall_2_1; }
		
		//properties+=Property*
		public Assignment getPropertiesAssignment_2_2() { return cPropertiesAssignment_2_2; }
		
		//Property
		public RuleCall getPropertiesPropertyParserRuleCall_2_2_0() { return cPropertiesPropertyParserRuleCall_2_2_0; }
		
		//labels+=ElkLabel*
		public Assignment getLabelsAssignment_2_3() { return cLabelsAssignment_2_3; }
		
		//ElkLabel
		public RuleCall getLabelsElkLabelParserRuleCall_2_3_0() { return cLabelsElkLabelParserRuleCall_2_3_0; }
		
		//'}'
		public Keyword getRightCurlyBracketKeyword_2_4() { return cRightCurlyBracketKeyword_2_4; }
	}
	public class ShapeLayoutElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.text.ElkGraph.ShapeLayout");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Keyword cLayoutKeyword_0 = (Keyword)cGroup.eContents().get(0);
		private final Keyword cLeftSquareBracketKeyword_1 = (Keyword)cGroup.eContents().get(1);
		private final UnorderedGroup cUnorderedGroup_2 = (UnorderedGroup)cGroup.eContents().get(2);
		private final Group cGroup_2_0 = (Group)cUnorderedGroup_2.eContents().get(0);
		private final Keyword cPositionKeyword_2_0_0 = (Keyword)cGroup_2_0.eContents().get(0);
		private final Keyword cColonKeyword_2_0_1 = (Keyword)cGroup_2_0.eContents().get(1);
		private final Assignment cXAssignment_2_0_2 = (Assignment)cGroup_2_0.eContents().get(2);
		private final RuleCall cXNumberParserRuleCall_2_0_2_0 = (RuleCall)cXAssignment_2_0_2.eContents().get(0);
		private final Keyword cCommaKeyword_2_0_3 = (Keyword)cGroup_2_0.eContents().get(3);
		private final Assignment cYAssignment_2_0_4 = (Assignment)cGroup_2_0.eContents().get(4);
		private final RuleCall cYNumberParserRuleCall_2_0_4_0 = (RuleCall)cYAssignment_2_0_4.eContents().get(0);
		private final Group cGroup_2_1 = (Group)cUnorderedGroup_2.eContents().get(1);
		private final Keyword cSizeKeyword_2_1_0 = (Keyword)cGroup_2_1.eContents().get(0);
		private final Keyword cColonKeyword_2_1_1 = (Keyword)cGroup_2_1.eContents().get(1);
		private final Assignment cWidthAssignment_2_1_2 = (Assignment)cGroup_2_1.eContents().get(2);
		private final RuleCall cWidthNumberParserRuleCall_2_1_2_0 = (RuleCall)cWidthAssignment_2_1_2.eContents().get(0);
		private final Keyword cCommaKeyword_2_1_3 = (Keyword)cGroup_2_1.eContents().get(3);
		private final Assignment cHeightAssignment_2_1_4 = (Assignment)cGroup_2_1.eContents().get(4);
		private final RuleCall cHeightNumberParserRuleCall_2_1_4_0 = (RuleCall)cHeightAssignment_2_1_4.eContents().get(0);
		private final Keyword cRightSquareBracketKeyword_3 = (Keyword)cGroup.eContents().get(3);
		
		//fragment ShapeLayout returns ElkShape:
		//	'layout' '[' (('position' ':' x=Number ',' y=Number)?
		//	& ('size' ':' width=Number ',' height=Number)?) ']';
		@Override public ParserRule getRule() { return rule; }
		
		//'layout' '[' (('position' ':' x=Number ',' y=Number)? & ('size' ':' width=Number ',' height=Number)?) ']'
		public Group getGroup() { return cGroup; }
		
		//'layout'
		public Keyword getLayoutKeyword_0() { return cLayoutKeyword_0; }
		
		//'['
		public Keyword getLeftSquareBracketKeyword_1() { return cLeftSquareBracketKeyword_1; }
		
		//('position' ':' x=Number ',' y=Number)? & ('size' ':' width=Number ',' height=Number)?
		public UnorderedGroup getUnorderedGroup_2() { return cUnorderedGroup_2; }
		
		//('position' ':' x=Number ',' y=Number)?
		public Group getGroup_2_0() { return cGroup_2_0; }
		
		//'position'
		public Keyword getPositionKeyword_2_0_0() { return cPositionKeyword_2_0_0; }
		
		//':'
		public Keyword getColonKeyword_2_0_1() { return cColonKeyword_2_0_1; }
		
		//x=Number
		public Assignment getXAssignment_2_0_2() { return cXAssignment_2_0_2; }
		
		//Number
		public RuleCall getXNumberParserRuleCall_2_0_2_0() { return cXNumberParserRuleCall_2_0_2_0; }
		
		//','
		public Keyword getCommaKeyword_2_0_3() { return cCommaKeyword_2_0_3; }
		
		//y=Number
		public Assignment getYAssignment_2_0_4() { return cYAssignment_2_0_4; }
		
		//Number
		public RuleCall getYNumberParserRuleCall_2_0_4_0() { return cYNumberParserRuleCall_2_0_4_0; }
		
		//('size' ':' width=Number ',' height=Number)?
		public Group getGroup_2_1() { return cGroup_2_1; }
		
		//'size'
		public Keyword getSizeKeyword_2_1_0() { return cSizeKeyword_2_1_0; }
		
		//':'
		public Keyword getColonKeyword_2_1_1() { return cColonKeyword_2_1_1; }
		
		//width=Number
		public Assignment getWidthAssignment_2_1_2() { return cWidthAssignment_2_1_2; }
		
		//Number
		public RuleCall getWidthNumberParserRuleCall_2_1_2_0() { return cWidthNumberParserRuleCall_2_1_2_0; }
		
		//','
		public Keyword getCommaKeyword_2_1_3() { return cCommaKeyword_2_1_3; }
		
		//height=Number
		public Assignment getHeightAssignment_2_1_4() { return cHeightAssignment_2_1_4; }
		
		//Number
		public RuleCall getHeightNumberParserRuleCall_2_1_4_0() { return cHeightNumberParserRuleCall_2_1_4_0; }
		
		//']'
		public Keyword getRightSquareBracketKeyword_3() { return cRightSquareBracketKeyword_3; }
	}
	public class ElkEdgeElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.text.ElkGraph.ElkEdge");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Keyword cEdgeKeyword_0 = (Keyword)cGroup.eContents().get(0);
		private final Group cGroup_1 = (Group)cGroup.eContents().get(1);
		private final Assignment cIdentifierAssignment_1_0 = (Assignment)cGroup_1.eContents().get(0);
		private final RuleCall cIdentifierIDTerminalRuleCall_1_0_0 = (RuleCall)cIdentifierAssignment_1_0.eContents().get(0);
		private final Keyword cColonKeyword_1_1 = (Keyword)cGroup_1.eContents().get(1);
		private final Assignment cSourcesAssignment_2 = (Assignment)cGroup.eContents().get(2);
		private final CrossReference cSourcesElkConnectableShapeCrossReference_2_0 = (CrossReference)cSourcesAssignment_2.eContents().get(0);
		private final RuleCall cSourcesElkConnectableShapeQualifiedIdParserRuleCall_2_0_1 = (RuleCall)cSourcesElkConnectableShapeCrossReference_2_0.eContents().get(1);
		private final Group cGroup_3 = (Group)cGroup.eContents().get(3);
		private final Keyword cCommaKeyword_3_0 = (Keyword)cGroup_3.eContents().get(0);
		private final Assignment cSourcesAssignment_3_1 = (Assignment)cGroup_3.eContents().get(1);
		private final CrossReference cSourcesElkConnectableShapeCrossReference_3_1_0 = (CrossReference)cSourcesAssignment_3_1.eContents().get(0);
		private final RuleCall cSourcesElkConnectableShapeQualifiedIdParserRuleCall_3_1_0_1 = (RuleCall)cSourcesElkConnectableShapeCrossReference_3_1_0.eContents().get(1);
		private final Keyword cHyphenMinusGreaterThanSignKeyword_4 = (Keyword)cGroup.eContents().get(4);
		private final Assignment cTargetsAssignment_5 = (Assignment)cGroup.eContents().get(5);
		private final CrossReference cTargetsElkConnectableShapeCrossReference_5_0 = (CrossReference)cTargetsAssignment_5.eContents().get(0);
		private final RuleCall cTargetsElkConnectableShapeQualifiedIdParserRuleCall_5_0_1 = (RuleCall)cTargetsElkConnectableShapeCrossReference_5_0.eContents().get(1);
		private final Group cGroup_6 = (Group)cGroup.eContents().get(6);
		private final Keyword cCommaKeyword_6_0 = (Keyword)cGroup_6.eContents().get(0);
		private final Assignment cTargetsAssignment_6_1 = (Assignment)cGroup_6.eContents().get(1);
		private final CrossReference cTargetsElkConnectableShapeCrossReference_6_1_0 = (CrossReference)cTargetsAssignment_6_1.eContents().get(0);
		private final RuleCall cTargetsElkConnectableShapeQualifiedIdParserRuleCall_6_1_0_1 = (RuleCall)cTargetsElkConnectableShapeCrossReference_6_1_0.eContents().get(1);
		private final Group cGroup_7 = (Group)cGroup.eContents().get(7);
		private final Keyword cLeftCurlyBracketKeyword_7_0 = (Keyword)cGroup_7.eContents().get(0);
		private final RuleCall cEdgeLayoutParserRuleCall_7_1 = (RuleCall)cGroup_7.eContents().get(1);
		private final Assignment cPropertiesAssignment_7_2 = (Assignment)cGroup_7.eContents().get(2);
		private final RuleCall cPropertiesPropertyParserRuleCall_7_2_0 = (RuleCall)cPropertiesAssignment_7_2.eContents().get(0);
		private final Assignment cLabelsAssignment_7_3 = (Assignment)cGroup_7.eContents().get(3);
		private final RuleCall cLabelsElkLabelParserRuleCall_7_3_0 = (RuleCall)cLabelsAssignment_7_3.eContents().get(0);
		private final Keyword cRightCurlyBracketKeyword_7_4 = (Keyword)cGroup_7.eContents().get(4);
		
		///* SuppressWarnings[BidirectionalReference] */ ElkEdge:
		//	'edge' (identifier=ID ':')?
		//	sources+=[ElkConnectableShape|QualifiedId] (',' sources+=[ElkConnectableShape|QualifiedId])* '->'
		//	targets+=[ElkConnectableShape|QualifiedId] (',' targets+=[ElkConnectableShape|QualifiedId])* ('{'
		//	EdgeLayout?
		//	properties+=Property*
		//	labels+=ElkLabel*
		//	'}')?;
		@Override public ParserRule getRule() { return rule; }
		
		//'edge' (identifier=ID ':')? sources+=[ElkConnectableShape|QualifiedId] (',' sources+=[ElkConnectableShape|QualifiedId])*
		//'->' targets+=[ElkConnectableShape|QualifiedId] (',' targets+=[ElkConnectableShape|QualifiedId])* ('{' EdgeLayout?
		//properties+=Property* labels+=ElkLabel* '}')?
		public Group getGroup() { return cGroup; }
		
		//'edge'
		public Keyword getEdgeKeyword_0() { return cEdgeKeyword_0; }
		
		//(identifier=ID ':')?
		public Group getGroup_1() { return cGroup_1; }
		
		//identifier=ID
		public Assignment getIdentifierAssignment_1_0() { return cIdentifierAssignment_1_0; }
		
		//ID
		public RuleCall getIdentifierIDTerminalRuleCall_1_0_0() { return cIdentifierIDTerminalRuleCall_1_0_0; }
		
		//':'
		public Keyword getColonKeyword_1_1() { return cColonKeyword_1_1; }
		
		//sources+=[ElkConnectableShape|QualifiedId]
		public Assignment getSourcesAssignment_2() { return cSourcesAssignment_2; }
		
		//[ElkConnectableShape|QualifiedId]
		public CrossReference getSourcesElkConnectableShapeCrossReference_2_0() { return cSourcesElkConnectableShapeCrossReference_2_0; }
		
		//QualifiedId
		public RuleCall getSourcesElkConnectableShapeQualifiedIdParserRuleCall_2_0_1() { return cSourcesElkConnectableShapeQualifiedIdParserRuleCall_2_0_1; }
		
		//(',' sources+=[ElkConnectableShape|QualifiedId])*
		public Group getGroup_3() { return cGroup_3; }
		
		//','
		public Keyword getCommaKeyword_3_0() { return cCommaKeyword_3_0; }
		
		//sources+=[ElkConnectableShape|QualifiedId]
		public Assignment getSourcesAssignment_3_1() { return cSourcesAssignment_3_1; }
		
		//[ElkConnectableShape|QualifiedId]
		public CrossReference getSourcesElkConnectableShapeCrossReference_3_1_0() { return cSourcesElkConnectableShapeCrossReference_3_1_0; }
		
		//QualifiedId
		public RuleCall getSourcesElkConnectableShapeQualifiedIdParserRuleCall_3_1_0_1() { return cSourcesElkConnectableShapeQualifiedIdParserRuleCall_3_1_0_1; }
		
		//'->'
		public Keyword getHyphenMinusGreaterThanSignKeyword_4() { return cHyphenMinusGreaterThanSignKeyword_4; }
		
		//targets+=[ElkConnectableShape|QualifiedId]
		public Assignment getTargetsAssignment_5() { return cTargetsAssignment_5; }
		
		//[ElkConnectableShape|QualifiedId]
		public CrossReference getTargetsElkConnectableShapeCrossReference_5_0() { return cTargetsElkConnectableShapeCrossReference_5_0; }
		
		//QualifiedId
		public RuleCall getTargetsElkConnectableShapeQualifiedIdParserRuleCall_5_0_1() { return cTargetsElkConnectableShapeQualifiedIdParserRuleCall_5_0_1; }
		
		//(',' targets+=[ElkConnectableShape|QualifiedId])*
		public Group getGroup_6() { return cGroup_6; }
		
		//','
		public Keyword getCommaKeyword_6_0() { return cCommaKeyword_6_0; }
		
		//targets+=[ElkConnectableShape|QualifiedId]
		public Assignment getTargetsAssignment_6_1() { return cTargetsAssignment_6_1; }
		
		//[ElkConnectableShape|QualifiedId]
		public CrossReference getTargetsElkConnectableShapeCrossReference_6_1_0() { return cTargetsElkConnectableShapeCrossReference_6_1_0; }
		
		//QualifiedId
		public RuleCall getTargetsElkConnectableShapeQualifiedIdParserRuleCall_6_1_0_1() { return cTargetsElkConnectableShapeQualifiedIdParserRuleCall_6_1_0_1; }
		
		//('{' EdgeLayout? properties+=Property* labels+=ElkLabel* '}')?
		public Group getGroup_7() { return cGroup_7; }
		
		//'{'
		public Keyword getLeftCurlyBracketKeyword_7_0() { return cLeftCurlyBracketKeyword_7_0; }
		
		//EdgeLayout?
		public RuleCall getEdgeLayoutParserRuleCall_7_1() { return cEdgeLayoutParserRuleCall_7_1; }
		
		//properties+=Property*
		public Assignment getPropertiesAssignment_7_2() { return cPropertiesAssignment_7_2; }
		
		//Property
		public RuleCall getPropertiesPropertyParserRuleCall_7_2_0() { return cPropertiesPropertyParserRuleCall_7_2_0; }
		
		//labels+=ElkLabel*
		public Assignment getLabelsAssignment_7_3() { return cLabelsAssignment_7_3; }
		
		//ElkLabel
		public RuleCall getLabelsElkLabelParserRuleCall_7_3_0() { return cLabelsElkLabelParserRuleCall_7_3_0; }
		
		//'}'
		public Keyword getRightCurlyBracketKeyword_7_4() { return cRightCurlyBracketKeyword_7_4; }
	}
	public class EdgeLayoutElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.text.ElkGraph.EdgeLayout");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Keyword cLayoutKeyword_0 = (Keyword)cGroup.eContents().get(0);
		private final Keyword cLeftSquareBracketKeyword_1 = (Keyword)cGroup.eContents().get(1);
		private final Alternatives cAlternatives_2 = (Alternatives)cGroup.eContents().get(2);
		private final Assignment cSectionsAssignment_2_0 = (Assignment)cAlternatives_2.eContents().get(0);
		private final RuleCall cSectionsElkSingleEdgeSectionParserRuleCall_2_0_0 = (RuleCall)cSectionsAssignment_2_0.eContents().get(0);
		private final Assignment cSectionsAssignment_2_1 = (Assignment)cAlternatives_2.eContents().get(1);
		private final RuleCall cSectionsElkEdgeSectionParserRuleCall_2_1_0 = (RuleCall)cSectionsAssignment_2_1.eContents().get(0);
		private final Keyword cRightSquareBracketKeyword_3 = (Keyword)cGroup.eContents().get(3);
		
		//fragment EdgeLayout returns ElkEdge:
		//	'layout' '[' (sections+=ElkSingleEdgeSection | sections+=ElkEdgeSection+) ']';
		@Override public ParserRule getRule() { return rule; }
		
		//'layout' '[' (sections+=ElkSingleEdgeSection | sections+=ElkEdgeSection+) ']'
		public Group getGroup() { return cGroup; }
		
		//'layout'
		public Keyword getLayoutKeyword_0() { return cLayoutKeyword_0; }
		
		//'['
		public Keyword getLeftSquareBracketKeyword_1() { return cLeftSquareBracketKeyword_1; }
		
		//sections+=ElkSingleEdgeSection | sections+=ElkEdgeSection+
		public Alternatives getAlternatives_2() { return cAlternatives_2; }
		
		//sections+=ElkSingleEdgeSection
		public Assignment getSectionsAssignment_2_0() { return cSectionsAssignment_2_0; }
		
		//ElkSingleEdgeSection
		public RuleCall getSectionsElkSingleEdgeSectionParserRuleCall_2_0_0() { return cSectionsElkSingleEdgeSectionParserRuleCall_2_0_0; }
		
		//sections+=ElkEdgeSection+
		public Assignment getSectionsAssignment_2_1() { return cSectionsAssignment_2_1; }
		
		//ElkEdgeSection
		public RuleCall getSectionsElkEdgeSectionParserRuleCall_2_1_0() { return cSectionsElkEdgeSectionParserRuleCall_2_1_0; }
		
		//']'
		public Keyword getRightSquareBracketKeyword_3() { return cRightSquareBracketKeyword_3; }
	}
	public class ElkSingleEdgeSectionElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.text.ElkGraph.ElkSingleEdgeSection");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Action cElkEdgeSectionAction_0 = (Action)cGroup.eContents().get(0);
		private final Group cGroup_1 = (Group)cGroup.eContents().get(1);
		private final UnorderedGroup cUnorderedGroup_1_0 = (UnorderedGroup)cGroup_1.eContents().get(0);
		private final Group cGroup_1_0_0 = (Group)cUnorderedGroup_1_0.eContents().get(0);
		private final Keyword cIncomingKeyword_1_0_0_0 = (Keyword)cGroup_1_0_0.eContents().get(0);
		private final Keyword cColonKeyword_1_0_0_1 = (Keyword)cGroup_1_0_0.eContents().get(1);
		private final Assignment cIncomingShapeAssignment_1_0_0_2 = (Assignment)cGroup_1_0_0.eContents().get(2);
		private final CrossReference cIncomingShapeElkConnectableShapeCrossReference_1_0_0_2_0 = (CrossReference)cIncomingShapeAssignment_1_0_0_2.eContents().get(0);
		private final RuleCall cIncomingShapeElkConnectableShapeQualifiedIdParserRuleCall_1_0_0_2_0_1 = (RuleCall)cIncomingShapeElkConnectableShapeCrossReference_1_0_0_2_0.eContents().get(1);
		private final Group cGroup_1_0_1 = (Group)cUnorderedGroup_1_0.eContents().get(1);
		private final Keyword cOutgoingKeyword_1_0_1_0 = (Keyword)cGroup_1_0_1.eContents().get(0);
		private final Keyword cColonKeyword_1_0_1_1 = (Keyword)cGroup_1_0_1.eContents().get(1);
		private final Assignment cOutgoingShapeAssignment_1_0_1_2 = (Assignment)cGroup_1_0_1.eContents().get(2);
		private final CrossReference cOutgoingShapeElkConnectableShapeCrossReference_1_0_1_2_0 = (CrossReference)cOutgoingShapeAssignment_1_0_1_2.eContents().get(0);
		private final RuleCall cOutgoingShapeElkConnectableShapeQualifiedIdParserRuleCall_1_0_1_2_0_1 = (RuleCall)cOutgoingShapeElkConnectableShapeCrossReference_1_0_1_2_0.eContents().get(1);
		private final Group cGroup_1_0_2 = (Group)cUnorderedGroup_1_0.eContents().get(2);
		private final Keyword cStartKeyword_1_0_2_0 = (Keyword)cGroup_1_0_2.eContents().get(0);
		private final Keyword cColonKeyword_1_0_2_1 = (Keyword)cGroup_1_0_2.eContents().get(1);
		private final Assignment cStartXAssignment_1_0_2_2 = (Assignment)cGroup_1_0_2.eContents().get(2);
		private final RuleCall cStartXNumberParserRuleCall_1_0_2_2_0 = (RuleCall)cStartXAssignment_1_0_2_2.eContents().get(0);
		private final Keyword cCommaKeyword_1_0_2_3 = (Keyword)cGroup_1_0_2.eContents().get(3);
		private final Assignment cStartYAssignment_1_0_2_4 = (Assignment)cGroup_1_0_2.eContents().get(4);
		private final RuleCall cStartYNumberParserRuleCall_1_0_2_4_0 = (RuleCall)cStartYAssignment_1_0_2_4.eContents().get(0);
		private final Group cGroup_1_0_3 = (Group)cUnorderedGroup_1_0.eContents().get(3);
		private final Keyword cEndKeyword_1_0_3_0 = (Keyword)cGroup_1_0_3.eContents().get(0);
		private final Keyword cColonKeyword_1_0_3_1 = (Keyword)cGroup_1_0_3.eContents().get(1);
		private final Assignment cEndXAssignment_1_0_3_2 = (Assignment)cGroup_1_0_3.eContents().get(2);
		private final RuleCall cEndXNumberParserRuleCall_1_0_3_2_0 = (RuleCall)cEndXAssignment_1_0_3_2.eContents().get(0);
		private final Keyword cCommaKeyword_1_0_3_3 = (Keyword)cGroup_1_0_3.eContents().get(3);
		private final Assignment cEndYAssignment_1_0_3_4 = (Assignment)cGroup_1_0_3.eContents().get(4);
		private final RuleCall cEndYNumberParserRuleCall_1_0_3_4_0 = (RuleCall)cEndYAssignment_1_0_3_4.eContents().get(0);
		private final Group cGroup_1_1 = (Group)cGroup_1.eContents().get(1);
		private final Keyword cBendsKeyword_1_1_0 = (Keyword)cGroup_1_1.eContents().get(0);
		private final Keyword cColonKeyword_1_1_1 = (Keyword)cGroup_1_1.eContents().get(1);
		private final Assignment cBendPointsAssignment_1_1_2 = (Assignment)cGroup_1_1.eContents().get(2);
		private final RuleCall cBendPointsElkBendPointParserRuleCall_1_1_2_0 = (RuleCall)cBendPointsAssignment_1_1_2.eContents().get(0);
		private final Group cGroup_1_1_3 = (Group)cGroup_1_1.eContents().get(3);
		private final Keyword cVerticalLineKeyword_1_1_3_0 = (Keyword)cGroup_1_1_3.eContents().get(0);
		private final Assignment cBendPointsAssignment_1_1_3_1 = (Assignment)cGroup_1_1_3.eContents().get(1);
		private final RuleCall cBendPointsElkBendPointParserRuleCall_1_1_3_1_0 = (RuleCall)cBendPointsAssignment_1_1_3_1.eContents().get(0);
		private final Assignment cPropertiesAssignment_1_2 = (Assignment)cGroup_1.eContents().get(2);
		private final RuleCall cPropertiesPropertyParserRuleCall_1_2_0 = (RuleCall)cPropertiesAssignment_1_2.eContents().get(0);
		
		//ElkSingleEdgeSection ElkEdgeSection:
		//	{ElkEdgeSection} ((('incoming' ':' incomingShape=[ElkConnectableShape|QualifiedId])?
		//	& ('outgoing' ':' outgoingShape=[ElkConnectableShape|QualifiedId])?
		//	& ('start' ':' startX=Number ',' startY=Number)?
		//	& ('end' ':' endX=Number ',' endY=Number)?) ('bends' ':' bendPoints+=ElkBendPoint ('|' bendPoints+=ElkBendPoint)*)?
		//	properties+=Property*);
		@Override public ParserRule getRule() { return rule; }
		
		//{ElkEdgeSection} ((('incoming' ':' incomingShape=[ElkConnectableShape|QualifiedId])? & ('outgoing' ':'
		//outgoingShape=[ElkConnectableShape|QualifiedId])? & ('start' ':' startX=Number ',' startY=Number)? & ('end' ':'
		//endX=Number ',' endY=Number)?) ('bends' ':' bendPoints+=ElkBendPoint ('|' bendPoints+=ElkBendPoint)*)?
		//properties+=Property*)
		public Group getGroup() { return cGroup; }
		
		//{ElkEdgeSection}
		public Action getElkEdgeSectionAction_0() { return cElkEdgeSectionAction_0; }
		
		//(('incoming' ':' incomingShape=[ElkConnectableShape|QualifiedId])? & ('outgoing' ':'
		//outgoingShape=[ElkConnectableShape|QualifiedId])? & ('start' ':' startX=Number ',' startY=Number)? & ('end' ':'
		//endX=Number ',' endY=Number)?) ('bends' ':' bendPoints+=ElkBendPoint ('|' bendPoints+=ElkBendPoint)*)?
		//properties+=Property*
		public Group getGroup_1() { return cGroup_1; }
		
		//('incoming' ':' incomingShape=[ElkConnectableShape|QualifiedId])? & ('outgoing' ':'
		//outgoingShape=[ElkConnectableShape|QualifiedId])? & ('start' ':' startX=Number ',' startY=Number)? & ('end' ':'
		//endX=Number ',' endY=Number)?
		public UnorderedGroup getUnorderedGroup_1_0() { return cUnorderedGroup_1_0; }
		
		//('incoming' ':' incomingShape=[ElkConnectableShape|QualifiedId])?
		public Group getGroup_1_0_0() { return cGroup_1_0_0; }
		
		//'incoming'
		public Keyword getIncomingKeyword_1_0_0_0() { return cIncomingKeyword_1_0_0_0; }
		
		//':'
		public Keyword getColonKeyword_1_0_0_1() { return cColonKeyword_1_0_0_1; }
		
		//incomingShape=[ElkConnectableShape|QualifiedId]
		public Assignment getIncomingShapeAssignment_1_0_0_2() { return cIncomingShapeAssignment_1_0_0_2; }
		
		//[ElkConnectableShape|QualifiedId]
		public CrossReference getIncomingShapeElkConnectableShapeCrossReference_1_0_0_2_0() { return cIncomingShapeElkConnectableShapeCrossReference_1_0_0_2_0; }
		
		//QualifiedId
		public RuleCall getIncomingShapeElkConnectableShapeQualifiedIdParserRuleCall_1_0_0_2_0_1() { return cIncomingShapeElkConnectableShapeQualifiedIdParserRuleCall_1_0_0_2_0_1; }
		
		//('outgoing' ':' outgoingShape=[ElkConnectableShape|QualifiedId])?
		public Group getGroup_1_0_1() { return cGroup_1_0_1; }
		
		//'outgoing'
		public Keyword getOutgoingKeyword_1_0_1_0() { return cOutgoingKeyword_1_0_1_0; }
		
		//':'
		public Keyword getColonKeyword_1_0_1_1() { return cColonKeyword_1_0_1_1; }
		
		//outgoingShape=[ElkConnectableShape|QualifiedId]
		public Assignment getOutgoingShapeAssignment_1_0_1_2() { return cOutgoingShapeAssignment_1_0_1_2; }
		
		//[ElkConnectableShape|QualifiedId]
		public CrossReference getOutgoingShapeElkConnectableShapeCrossReference_1_0_1_2_0() { return cOutgoingShapeElkConnectableShapeCrossReference_1_0_1_2_0; }
		
		//QualifiedId
		public RuleCall getOutgoingShapeElkConnectableShapeQualifiedIdParserRuleCall_1_0_1_2_0_1() { return cOutgoingShapeElkConnectableShapeQualifiedIdParserRuleCall_1_0_1_2_0_1; }
		
		//('start' ':' startX=Number ',' startY=Number)?
		public Group getGroup_1_0_2() { return cGroup_1_0_2; }
		
		//'start'
		public Keyword getStartKeyword_1_0_2_0() { return cStartKeyword_1_0_2_0; }
		
		//':'
		public Keyword getColonKeyword_1_0_2_1() { return cColonKeyword_1_0_2_1; }
		
		//startX=Number
		public Assignment getStartXAssignment_1_0_2_2() { return cStartXAssignment_1_0_2_2; }
		
		//Number
		public RuleCall getStartXNumberParserRuleCall_1_0_2_2_0() { return cStartXNumberParserRuleCall_1_0_2_2_0; }
		
		//','
		public Keyword getCommaKeyword_1_0_2_3() { return cCommaKeyword_1_0_2_3; }
		
		//startY=Number
		public Assignment getStartYAssignment_1_0_2_4() { return cStartYAssignment_1_0_2_4; }
		
		//Number
		public RuleCall getStartYNumberParserRuleCall_1_0_2_4_0() { return cStartYNumberParserRuleCall_1_0_2_4_0; }
		
		//('end' ':' endX=Number ',' endY=Number)?
		public Group getGroup_1_0_3() { return cGroup_1_0_3; }
		
		//'end'
		public Keyword getEndKeyword_1_0_3_0() { return cEndKeyword_1_0_3_0; }
		
		//':'
		public Keyword getColonKeyword_1_0_3_1() { return cColonKeyword_1_0_3_1; }
		
		//endX=Number
		public Assignment getEndXAssignment_1_0_3_2() { return cEndXAssignment_1_0_3_2; }
		
		//Number
		public RuleCall getEndXNumberParserRuleCall_1_0_3_2_0() { return cEndXNumberParserRuleCall_1_0_3_2_0; }
		
		//','
		public Keyword getCommaKeyword_1_0_3_3() { return cCommaKeyword_1_0_3_3; }
		
		//endY=Number
		public Assignment getEndYAssignment_1_0_3_4() { return cEndYAssignment_1_0_3_4; }
		
		//Number
		public RuleCall getEndYNumberParserRuleCall_1_0_3_4_0() { return cEndYNumberParserRuleCall_1_0_3_4_0; }
		
		//('bends' ':' bendPoints+=ElkBendPoint ('|' bendPoints+=ElkBendPoint)*)?
		public Group getGroup_1_1() { return cGroup_1_1; }
		
		//'bends'
		public Keyword getBendsKeyword_1_1_0() { return cBendsKeyword_1_1_0; }
		
		//':'
		public Keyword getColonKeyword_1_1_1() { return cColonKeyword_1_1_1; }
		
		//bendPoints+=ElkBendPoint
		public Assignment getBendPointsAssignment_1_1_2() { return cBendPointsAssignment_1_1_2; }
		
		//ElkBendPoint
		public RuleCall getBendPointsElkBendPointParserRuleCall_1_1_2_0() { return cBendPointsElkBendPointParserRuleCall_1_1_2_0; }
		
		//('|' bendPoints+=ElkBendPoint)*
		public Group getGroup_1_1_3() { return cGroup_1_1_3; }
		
		//'|'
		public Keyword getVerticalLineKeyword_1_1_3_0() { return cVerticalLineKeyword_1_1_3_0; }
		
		//bendPoints+=ElkBendPoint
		public Assignment getBendPointsAssignment_1_1_3_1() { return cBendPointsAssignment_1_1_3_1; }
		
		//ElkBendPoint
		public RuleCall getBendPointsElkBendPointParserRuleCall_1_1_3_1_0() { return cBendPointsElkBendPointParserRuleCall_1_1_3_1_0; }
		
		//properties+=Property*
		public Assignment getPropertiesAssignment_1_2() { return cPropertiesAssignment_1_2; }
		
		//Property
		public RuleCall getPropertiesPropertyParserRuleCall_1_2_0() { return cPropertiesPropertyParserRuleCall_1_2_0; }
	}
	public class ElkEdgeSectionElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.text.ElkGraph.ElkEdgeSection");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Keyword cSectionKeyword_0 = (Keyword)cGroup.eContents().get(0);
		private final Assignment cIdentifierAssignment_1 = (Assignment)cGroup.eContents().get(1);
		private final RuleCall cIdentifierIDTerminalRuleCall_1_0 = (RuleCall)cIdentifierAssignment_1.eContents().get(0);
		private final Group cGroup_2 = (Group)cGroup.eContents().get(2);
		private final Keyword cHyphenMinusGreaterThanSignKeyword_2_0 = (Keyword)cGroup_2.eContents().get(0);
		private final Assignment cOutgoingSectionsAssignment_2_1 = (Assignment)cGroup_2.eContents().get(1);
		private final CrossReference cOutgoingSectionsElkEdgeSectionCrossReference_2_1_0 = (CrossReference)cOutgoingSectionsAssignment_2_1.eContents().get(0);
		private final RuleCall cOutgoingSectionsElkEdgeSectionIDTerminalRuleCall_2_1_0_1 = (RuleCall)cOutgoingSectionsElkEdgeSectionCrossReference_2_1_0.eContents().get(1);
		private final Group cGroup_2_2 = (Group)cGroup_2.eContents().get(2);
		private final Keyword cCommaKeyword_2_2_0 = (Keyword)cGroup_2_2.eContents().get(0);
		private final Assignment cOutgoingSectionsAssignment_2_2_1 = (Assignment)cGroup_2_2.eContents().get(1);
		private final CrossReference cOutgoingSectionsElkEdgeSectionCrossReference_2_2_1_0 = (CrossReference)cOutgoingSectionsAssignment_2_2_1.eContents().get(0);
		private final RuleCall cOutgoingSectionsElkEdgeSectionIDTerminalRuleCall_2_2_1_0_1 = (RuleCall)cOutgoingSectionsElkEdgeSectionCrossReference_2_2_1_0.eContents().get(1);
		private final Keyword cLeftSquareBracketKeyword_3 = (Keyword)cGroup.eContents().get(3);
		private final Group cGroup_4 = (Group)cGroup.eContents().get(4);
		private final UnorderedGroup cUnorderedGroup_4_0 = (UnorderedGroup)cGroup_4.eContents().get(0);
		private final Group cGroup_4_0_0 = (Group)cUnorderedGroup_4_0.eContents().get(0);
		private final Keyword cIncomingKeyword_4_0_0_0 = (Keyword)cGroup_4_0_0.eContents().get(0);
		private final Keyword cColonKeyword_4_0_0_1 = (Keyword)cGroup_4_0_0.eContents().get(1);
		private final Assignment cIncomingShapeAssignment_4_0_0_2 = (Assignment)cGroup_4_0_0.eContents().get(2);
		private final CrossReference cIncomingShapeElkConnectableShapeCrossReference_4_0_0_2_0 = (CrossReference)cIncomingShapeAssignment_4_0_0_2.eContents().get(0);
		private final RuleCall cIncomingShapeElkConnectableShapeQualifiedIdParserRuleCall_4_0_0_2_0_1 = (RuleCall)cIncomingShapeElkConnectableShapeCrossReference_4_0_0_2_0.eContents().get(1);
		private final Group cGroup_4_0_1 = (Group)cUnorderedGroup_4_0.eContents().get(1);
		private final Keyword cOutgoingKeyword_4_0_1_0 = (Keyword)cGroup_4_0_1.eContents().get(0);
		private final Keyword cColonKeyword_4_0_1_1 = (Keyword)cGroup_4_0_1.eContents().get(1);
		private final Assignment cOutgoingShapeAssignment_4_0_1_2 = (Assignment)cGroup_4_0_1.eContents().get(2);
		private final CrossReference cOutgoingShapeElkConnectableShapeCrossReference_4_0_1_2_0 = (CrossReference)cOutgoingShapeAssignment_4_0_1_2.eContents().get(0);
		private final RuleCall cOutgoingShapeElkConnectableShapeQualifiedIdParserRuleCall_4_0_1_2_0_1 = (RuleCall)cOutgoingShapeElkConnectableShapeCrossReference_4_0_1_2_0.eContents().get(1);
		private final Group cGroup_4_0_2 = (Group)cUnorderedGroup_4_0.eContents().get(2);
		private final Keyword cStartKeyword_4_0_2_0 = (Keyword)cGroup_4_0_2.eContents().get(0);
		private final Keyword cColonKeyword_4_0_2_1 = (Keyword)cGroup_4_0_2.eContents().get(1);
		private final Assignment cStartXAssignment_4_0_2_2 = (Assignment)cGroup_4_0_2.eContents().get(2);
		private final RuleCall cStartXNumberParserRuleCall_4_0_2_2_0 = (RuleCall)cStartXAssignment_4_0_2_2.eContents().get(0);
		private final Keyword cCommaKeyword_4_0_2_3 = (Keyword)cGroup_4_0_2.eContents().get(3);
		private final Assignment cStartYAssignment_4_0_2_4 = (Assignment)cGroup_4_0_2.eContents().get(4);
		private final RuleCall cStartYNumberParserRuleCall_4_0_2_4_0 = (RuleCall)cStartYAssignment_4_0_2_4.eContents().get(0);
		private final Group cGroup_4_0_3 = (Group)cUnorderedGroup_4_0.eContents().get(3);
		private final Keyword cEndKeyword_4_0_3_0 = (Keyword)cGroup_4_0_3.eContents().get(0);
		private final Keyword cColonKeyword_4_0_3_1 = (Keyword)cGroup_4_0_3.eContents().get(1);
		private final Assignment cEndXAssignment_4_0_3_2 = (Assignment)cGroup_4_0_3.eContents().get(2);
		private final RuleCall cEndXNumberParserRuleCall_4_0_3_2_0 = (RuleCall)cEndXAssignment_4_0_3_2.eContents().get(0);
		private final Keyword cCommaKeyword_4_0_3_3 = (Keyword)cGroup_4_0_3.eContents().get(3);
		private final Assignment cEndYAssignment_4_0_3_4 = (Assignment)cGroup_4_0_3.eContents().get(4);
		private final RuleCall cEndYNumberParserRuleCall_4_0_3_4_0 = (RuleCall)cEndYAssignment_4_0_3_4.eContents().get(0);
		private final Group cGroup_4_1 = (Group)cGroup_4.eContents().get(1);
		private final Keyword cBendsKeyword_4_1_0 = (Keyword)cGroup_4_1.eContents().get(0);
		private final Keyword cColonKeyword_4_1_1 = (Keyword)cGroup_4_1.eContents().get(1);
		private final Assignment cBendPointsAssignment_4_1_2 = (Assignment)cGroup_4_1.eContents().get(2);
		private final RuleCall cBendPointsElkBendPointParserRuleCall_4_1_2_0 = (RuleCall)cBendPointsAssignment_4_1_2.eContents().get(0);
		private final Group cGroup_4_1_3 = (Group)cGroup_4_1.eContents().get(3);
		private final Keyword cVerticalLineKeyword_4_1_3_0 = (Keyword)cGroup_4_1_3.eContents().get(0);
		private final Assignment cBendPointsAssignment_4_1_3_1 = (Assignment)cGroup_4_1_3.eContents().get(1);
		private final RuleCall cBendPointsElkBendPointParserRuleCall_4_1_3_1_0 = (RuleCall)cBendPointsAssignment_4_1_3_1.eContents().get(0);
		private final Assignment cPropertiesAssignment_4_2 = (Assignment)cGroup_4.eContents().get(2);
		private final RuleCall cPropertiesPropertyParserRuleCall_4_2_0 = (RuleCall)cPropertiesAssignment_4_2.eContents().get(0);
		private final Keyword cRightSquareBracketKeyword_5 = (Keyword)cGroup.eContents().get(5);
		
		///* SuppressWarnings[BidirectionalReference] */ ElkEdgeSection:
		//	'section' identifier=ID ('->' outgoingSections+=[ElkEdgeSection] (',' outgoingSections+=[ElkEdgeSection])*)? '['
		//	((('incoming' ':' incomingShape=[ElkConnectableShape|QualifiedId])?
		//	& ('outgoing' ':' outgoingShape=[ElkConnectableShape|QualifiedId])?
		//	& ('start' ':' startX=Number ',' startY=Number)?
		//	& ('end' ':' endX=Number ',' endY=Number)?) ('bends' ':' bendPoints+=ElkBendPoint ('|' bendPoints+=ElkBendPoint)*)?
		//	properties+=Property*) ']';
		@Override public ParserRule getRule() { return rule; }
		
		//'section' identifier=ID ('->' outgoingSections+=[ElkEdgeSection] (',' outgoingSections+=[ElkEdgeSection])*)? '['
		//((('incoming' ':' incomingShape=[ElkConnectableShape|QualifiedId])? & ('outgoing' ':'
		//outgoingShape=[ElkConnectableShape|QualifiedId])? & ('start' ':' startX=Number ',' startY=Number)? & ('end' ':'
		//endX=Number ',' endY=Number)?) ('bends' ':' bendPoints+=ElkBendPoint ('|' bendPoints+=ElkBendPoint)*)?
		//properties+=Property*) ']'
		public Group getGroup() { return cGroup; }
		
		//'section'
		public Keyword getSectionKeyword_0() { return cSectionKeyword_0; }
		
		//identifier=ID
		public Assignment getIdentifierAssignment_1() { return cIdentifierAssignment_1; }
		
		//ID
		public RuleCall getIdentifierIDTerminalRuleCall_1_0() { return cIdentifierIDTerminalRuleCall_1_0; }
		
		//('->' outgoingSections+=[ElkEdgeSection] (',' outgoingSections+=[ElkEdgeSection])*)?
		public Group getGroup_2() { return cGroup_2; }
		
		//'->'
		public Keyword getHyphenMinusGreaterThanSignKeyword_2_0() { return cHyphenMinusGreaterThanSignKeyword_2_0; }
		
		//outgoingSections+=[ElkEdgeSection]
		public Assignment getOutgoingSectionsAssignment_2_1() { return cOutgoingSectionsAssignment_2_1; }
		
		//[ElkEdgeSection]
		public CrossReference getOutgoingSectionsElkEdgeSectionCrossReference_2_1_0() { return cOutgoingSectionsElkEdgeSectionCrossReference_2_1_0; }
		
		//ID
		public RuleCall getOutgoingSectionsElkEdgeSectionIDTerminalRuleCall_2_1_0_1() { return cOutgoingSectionsElkEdgeSectionIDTerminalRuleCall_2_1_0_1; }
		
		//(',' outgoingSections+=[ElkEdgeSection])*
		public Group getGroup_2_2() { return cGroup_2_2; }
		
		//','
		public Keyword getCommaKeyword_2_2_0() { return cCommaKeyword_2_2_0; }
		
		//outgoingSections+=[ElkEdgeSection]
		public Assignment getOutgoingSectionsAssignment_2_2_1() { return cOutgoingSectionsAssignment_2_2_1; }
		
		//[ElkEdgeSection]
		public CrossReference getOutgoingSectionsElkEdgeSectionCrossReference_2_2_1_0() { return cOutgoingSectionsElkEdgeSectionCrossReference_2_2_1_0; }
		
		//ID
		public RuleCall getOutgoingSectionsElkEdgeSectionIDTerminalRuleCall_2_2_1_0_1() { return cOutgoingSectionsElkEdgeSectionIDTerminalRuleCall_2_2_1_0_1; }
		
		//'['
		public Keyword getLeftSquareBracketKeyword_3() { return cLeftSquareBracketKeyword_3; }
		
		//(('incoming' ':' incomingShape=[ElkConnectableShape|QualifiedId])? & ('outgoing' ':'
		//outgoingShape=[ElkConnectableShape|QualifiedId])? & ('start' ':' startX=Number ',' startY=Number)? & ('end' ':'
		//endX=Number ',' endY=Number)?) ('bends' ':' bendPoints+=ElkBendPoint ('|' bendPoints+=ElkBendPoint)*)?
		//properties+=Property*
		public Group getGroup_4() { return cGroup_4; }
		
		//('incoming' ':' incomingShape=[ElkConnectableShape|QualifiedId])? & ('outgoing' ':'
		//outgoingShape=[ElkConnectableShape|QualifiedId])? & ('start' ':' startX=Number ',' startY=Number)? & ('end' ':'
		//endX=Number ',' endY=Number)?
		public UnorderedGroup getUnorderedGroup_4_0() { return cUnorderedGroup_4_0; }
		
		//('incoming' ':' incomingShape=[ElkConnectableShape|QualifiedId])?
		public Group getGroup_4_0_0() { return cGroup_4_0_0; }
		
		//'incoming'
		public Keyword getIncomingKeyword_4_0_0_0() { return cIncomingKeyword_4_0_0_0; }
		
		//':'
		public Keyword getColonKeyword_4_0_0_1() { return cColonKeyword_4_0_0_1; }
		
		//incomingShape=[ElkConnectableShape|QualifiedId]
		public Assignment getIncomingShapeAssignment_4_0_0_2() { return cIncomingShapeAssignment_4_0_0_2; }
		
		//[ElkConnectableShape|QualifiedId]
		public CrossReference getIncomingShapeElkConnectableShapeCrossReference_4_0_0_2_0() { return cIncomingShapeElkConnectableShapeCrossReference_4_0_0_2_0; }
		
		//QualifiedId
		public RuleCall getIncomingShapeElkConnectableShapeQualifiedIdParserRuleCall_4_0_0_2_0_1() { return cIncomingShapeElkConnectableShapeQualifiedIdParserRuleCall_4_0_0_2_0_1; }
		
		//('outgoing' ':' outgoingShape=[ElkConnectableShape|QualifiedId])?
		public Group getGroup_4_0_1() { return cGroup_4_0_1; }
		
		//'outgoing'
		public Keyword getOutgoingKeyword_4_0_1_0() { return cOutgoingKeyword_4_0_1_0; }
		
		//':'
		public Keyword getColonKeyword_4_0_1_1() { return cColonKeyword_4_0_1_1; }
		
		//outgoingShape=[ElkConnectableShape|QualifiedId]
		public Assignment getOutgoingShapeAssignment_4_0_1_2() { return cOutgoingShapeAssignment_4_0_1_2; }
		
		//[ElkConnectableShape|QualifiedId]
		public CrossReference getOutgoingShapeElkConnectableShapeCrossReference_4_0_1_2_0() { return cOutgoingShapeElkConnectableShapeCrossReference_4_0_1_2_0; }
		
		//QualifiedId
		public RuleCall getOutgoingShapeElkConnectableShapeQualifiedIdParserRuleCall_4_0_1_2_0_1() { return cOutgoingShapeElkConnectableShapeQualifiedIdParserRuleCall_4_0_1_2_0_1; }
		
		//('start' ':' startX=Number ',' startY=Number)?
		public Group getGroup_4_0_2() { return cGroup_4_0_2; }
		
		//'start'
		public Keyword getStartKeyword_4_0_2_0() { return cStartKeyword_4_0_2_0; }
		
		//':'
		public Keyword getColonKeyword_4_0_2_1() { return cColonKeyword_4_0_2_1; }
		
		//startX=Number
		public Assignment getStartXAssignment_4_0_2_2() { return cStartXAssignment_4_0_2_2; }
		
		//Number
		public RuleCall getStartXNumberParserRuleCall_4_0_2_2_0() { return cStartXNumberParserRuleCall_4_0_2_2_0; }
		
		//','
		public Keyword getCommaKeyword_4_0_2_3() { return cCommaKeyword_4_0_2_3; }
		
		//startY=Number
		public Assignment getStartYAssignment_4_0_2_4() { return cStartYAssignment_4_0_2_4; }
		
		//Number
		public RuleCall getStartYNumberParserRuleCall_4_0_2_4_0() { return cStartYNumberParserRuleCall_4_0_2_4_0; }
		
		//('end' ':' endX=Number ',' endY=Number)?
		public Group getGroup_4_0_3() { return cGroup_4_0_3; }
		
		//'end'
		public Keyword getEndKeyword_4_0_3_0() { return cEndKeyword_4_0_3_0; }
		
		//':'
		public Keyword getColonKeyword_4_0_3_1() { return cColonKeyword_4_0_3_1; }
		
		//endX=Number
		public Assignment getEndXAssignment_4_0_3_2() { return cEndXAssignment_4_0_3_2; }
		
		//Number
		public RuleCall getEndXNumberParserRuleCall_4_0_3_2_0() { return cEndXNumberParserRuleCall_4_0_3_2_0; }
		
		//','
		public Keyword getCommaKeyword_4_0_3_3() { return cCommaKeyword_4_0_3_3; }
		
		//endY=Number
		public Assignment getEndYAssignment_4_0_3_4() { return cEndYAssignment_4_0_3_4; }
		
		//Number
		public RuleCall getEndYNumberParserRuleCall_4_0_3_4_0() { return cEndYNumberParserRuleCall_4_0_3_4_0; }
		
		//('bends' ':' bendPoints+=ElkBendPoint ('|' bendPoints+=ElkBendPoint)*)?
		public Group getGroup_4_1() { return cGroup_4_1; }
		
		//'bends'
		public Keyword getBendsKeyword_4_1_0() { return cBendsKeyword_4_1_0; }
		
		//':'
		public Keyword getColonKeyword_4_1_1() { return cColonKeyword_4_1_1; }
		
		//bendPoints+=ElkBendPoint
		public Assignment getBendPointsAssignment_4_1_2() { return cBendPointsAssignment_4_1_2; }
		
		//ElkBendPoint
		public RuleCall getBendPointsElkBendPointParserRuleCall_4_1_2_0() { return cBendPointsElkBendPointParserRuleCall_4_1_2_0; }
		
		//('|' bendPoints+=ElkBendPoint)*
		public Group getGroup_4_1_3() { return cGroup_4_1_3; }
		
		//'|'
		public Keyword getVerticalLineKeyword_4_1_3_0() { return cVerticalLineKeyword_4_1_3_0; }
		
		//bendPoints+=ElkBendPoint
		public Assignment getBendPointsAssignment_4_1_3_1() { return cBendPointsAssignment_4_1_3_1; }
		
		//ElkBendPoint
		public RuleCall getBendPointsElkBendPointParserRuleCall_4_1_3_1_0() { return cBendPointsElkBendPointParserRuleCall_4_1_3_1_0; }
		
		//properties+=Property*
		public Assignment getPropertiesAssignment_4_2() { return cPropertiesAssignment_4_2; }
		
		//Property
		public RuleCall getPropertiesPropertyParserRuleCall_4_2_0() { return cPropertiesPropertyParserRuleCall_4_2_0; }
		
		//']'
		public Keyword getRightSquareBracketKeyword_5() { return cRightSquareBracketKeyword_5; }
	}
	public class ElkBendPointElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.text.ElkGraph.ElkBendPoint");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Assignment cXAssignment_0 = (Assignment)cGroup.eContents().get(0);
		private final RuleCall cXNumberParserRuleCall_0_0 = (RuleCall)cXAssignment_0.eContents().get(0);
		private final Keyword cCommaKeyword_1 = (Keyword)cGroup.eContents().get(1);
		private final Assignment cYAssignment_2 = (Assignment)cGroup.eContents().get(2);
		private final RuleCall cYNumberParserRuleCall_2_0 = (RuleCall)cYAssignment_2.eContents().get(0);
		
		//ElkBendPoint:
		//	x=Number ',' y=Number;
		@Override public ParserRule getRule() { return rule; }
		
		//x=Number ',' y=Number
		public Group getGroup() { return cGroup; }
		
		//x=Number
		public Assignment getXAssignment_0() { return cXAssignment_0; }
		
		//Number
		public RuleCall getXNumberParserRuleCall_0_0() { return cXNumberParserRuleCall_0_0; }
		
		//','
		public Keyword getCommaKeyword_1() { return cCommaKeyword_1; }
		
		//y=Number
		public Assignment getYAssignment_2() { return cYAssignment_2; }
		
		//Number
		public RuleCall getYNumberParserRuleCall_2_0() { return cYNumberParserRuleCall_2_0; }
	}
	public class QualifiedIdElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.text.ElkGraph.QualifiedId");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final RuleCall cIDTerminalRuleCall_0 = (RuleCall)cGroup.eContents().get(0);
		private final Group cGroup_1 = (Group)cGroup.eContents().get(1);
		private final Keyword cFullStopKeyword_1_0 = (Keyword)cGroup_1.eContents().get(0);
		private final RuleCall cIDTerminalRuleCall_1_1 = (RuleCall)cGroup_1.eContents().get(1);
		
		//QualifiedId:
		//	ID ('.' ID)*;
		@Override public ParserRule getRule() { return rule; }
		
		//ID ('.' ID)*
		public Group getGroup() { return cGroup; }
		
		//ID
		public RuleCall getIDTerminalRuleCall_0() { return cIDTerminalRuleCall_0; }
		
		//('.' ID)*
		public Group getGroup_1() { return cGroup_1; }
		
		//'.'
		public Keyword getFullStopKeyword_1_0() { return cFullStopKeyword_1_0; }
		
		//ID
		public RuleCall getIDTerminalRuleCall_1_1() { return cIDTerminalRuleCall_1_1; }
	}
	public class NumberElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.text.ElkGraph.Number");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final RuleCall cSIGNED_INTTerminalRuleCall_0 = (RuleCall)cAlternatives.eContents().get(0);
		private final RuleCall cFLOATTerminalRuleCall_1 = (RuleCall)cAlternatives.eContents().get(1);
		
		//Number ecore::EDouble:
		//	SIGNED_INT | FLOAT;
		@Override public ParserRule getRule() { return rule; }
		
		//SIGNED_INT | FLOAT
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//SIGNED_INT
		public RuleCall getSIGNED_INTTerminalRuleCall_0() { return cSIGNED_INTTerminalRuleCall_0; }
		
		//FLOAT
		public RuleCall getFLOATTerminalRuleCall_1() { return cFLOATTerminalRuleCall_1; }
	}
	public class PropertyElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.text.ElkGraph.Property");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Assignment cKeyAssignment_0 = (Assignment)cGroup.eContents().get(0);
		private final RuleCall cKeyPropertyKeyParserRuleCall_0_0 = (RuleCall)cKeyAssignment_0.eContents().get(0);
		private final Keyword cColonKeyword_1 = (Keyword)cGroup.eContents().get(1);
		private final Alternatives cAlternatives_2 = (Alternatives)cGroup.eContents().get(2);
		private final Assignment cValueAssignment_2_0 = (Assignment)cAlternatives_2.eContents().get(0);
		private final RuleCall cValueStringValueParserRuleCall_2_0_0 = (RuleCall)cValueAssignment_2_0.eContents().get(0);
		private final Assignment cValueAssignment_2_1 = (Assignment)cAlternatives_2.eContents().get(1);
		private final RuleCall cValueQualifiedIdValueParserRuleCall_2_1_0 = (RuleCall)cValueAssignment_2_1.eContents().get(0);
		private final Assignment cValueAssignment_2_2 = (Assignment)cAlternatives_2.eContents().get(2);
		private final RuleCall cValueNumberValueParserRuleCall_2_2_0 = (RuleCall)cValueAssignment_2_2.eContents().get(0);
		private final Assignment cValueAssignment_2_3 = (Assignment)cAlternatives_2.eContents().get(3);
		private final RuleCall cValueBooleanValueParserRuleCall_2_3_0 = (RuleCall)cValueAssignment_2_3.eContents().get(0);
		private final Keyword cNullKeyword_2_4 = (Keyword)cAlternatives_2.eContents().get(4);
		
		//Property ElkPropertyToValueMapEntry:
		//	key=PropertyKey ':' (value=StringValue | value=QualifiedIdValue | value=NumberValue | value=BooleanValue | 'null');
		@Override public ParserRule getRule() { return rule; }
		
		//key=PropertyKey ':' (value=StringValue | value=QualifiedIdValue | value=NumberValue | value=BooleanValue | 'null')
		public Group getGroup() { return cGroup; }
		
		//key=PropertyKey
		public Assignment getKeyAssignment_0() { return cKeyAssignment_0; }
		
		//PropertyKey
		public RuleCall getKeyPropertyKeyParserRuleCall_0_0() { return cKeyPropertyKeyParserRuleCall_0_0; }
		
		//':'
		public Keyword getColonKeyword_1() { return cColonKeyword_1; }
		
		//value=StringValue | value=QualifiedIdValue | value=NumberValue | value=BooleanValue | 'null'
		public Alternatives getAlternatives_2() { return cAlternatives_2; }
		
		//value=StringValue
		public Assignment getValueAssignment_2_0() { return cValueAssignment_2_0; }
		
		//StringValue
		public RuleCall getValueStringValueParserRuleCall_2_0_0() { return cValueStringValueParserRuleCall_2_0_0; }
		
		//value=QualifiedIdValue
		public Assignment getValueAssignment_2_1() { return cValueAssignment_2_1; }
		
		//QualifiedIdValue
		public RuleCall getValueQualifiedIdValueParserRuleCall_2_1_0() { return cValueQualifiedIdValueParserRuleCall_2_1_0; }
		
		//value=NumberValue
		public Assignment getValueAssignment_2_2() { return cValueAssignment_2_2; }
		
		//NumberValue
		public RuleCall getValueNumberValueParserRuleCall_2_2_0() { return cValueNumberValueParserRuleCall_2_2_0; }
		
		//value=BooleanValue
		public Assignment getValueAssignment_2_3() { return cValueAssignment_2_3; }
		
		//BooleanValue
		public RuleCall getValueBooleanValueParserRuleCall_2_3_0() { return cValueBooleanValueParserRuleCall_2_3_0; }
		
		//'null'
		public Keyword getNullKeyword_2_4() { return cNullKeyword_2_4; }
	}
	public class PropertyKeyElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.text.ElkGraph.PropertyKey");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final RuleCall cIDTerminalRuleCall_0 = (RuleCall)cGroup.eContents().get(0);
		private final Group cGroup_1 = (Group)cGroup.eContents().get(1);
		private final Keyword cFullStopKeyword_1_0 = (Keyword)cGroup_1.eContents().get(0);
		private final RuleCall cIDTerminalRuleCall_1_1 = (RuleCall)cGroup_1.eContents().get(1);
		
		//PropertyKey IProperty hidden():
		//	ID ('.' ID)*;
		@Override public ParserRule getRule() { return rule; }
		
		//ID ('.' ID)*
		public Group getGroup() { return cGroup; }
		
		//ID
		public RuleCall getIDTerminalRuleCall_0() { return cIDTerminalRuleCall_0; }
		
		//('.' ID)*
		public Group getGroup_1() { return cGroup_1; }
		
		//'.'
		public Keyword getFullStopKeyword_1_0() { return cFullStopKeyword_1_0; }
		
		//ID
		public RuleCall getIDTerminalRuleCall_1_1() { return cIDTerminalRuleCall_1_1; }
	}
	public class StringValueElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.text.ElkGraph.StringValue");
		private final RuleCall cSTRINGTerminalRuleCall = (RuleCall)rule.eContents().get(1);
		
		//StringValue ecore::EJavaObject:
		//	STRING;
		@Override public ParserRule getRule() { return rule; }
		
		//STRING
		public RuleCall getSTRINGTerminalRuleCall() { return cSTRINGTerminalRuleCall; }
	}
	public class QualifiedIdValueElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.text.ElkGraph.QualifiedIdValue");
		private final RuleCall cQualifiedIdParserRuleCall = (RuleCall)rule.eContents().get(1);
		
		//QualifiedIdValue ecore::EJavaObject:
		//	QualifiedId;
		@Override public ParserRule getRule() { return rule; }
		
		//QualifiedId
		public RuleCall getQualifiedIdParserRuleCall() { return cQualifiedIdParserRuleCall; }
	}
	public class NumberValueElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.text.ElkGraph.NumberValue");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final RuleCall cSIGNED_INTTerminalRuleCall_0 = (RuleCall)cAlternatives.eContents().get(0);
		private final RuleCall cFLOATTerminalRuleCall_1 = (RuleCall)cAlternatives.eContents().get(1);
		
		//NumberValue ecore::EJavaObject:
		//	SIGNED_INT | FLOAT;
		@Override public ParserRule getRule() { return rule; }
		
		//SIGNED_INT | FLOAT
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//SIGNED_INT
		public RuleCall getSIGNED_INTTerminalRuleCall_0() { return cSIGNED_INTTerminalRuleCall_0; }
		
		//FLOAT
		public RuleCall getFLOATTerminalRuleCall_1() { return cFLOATTerminalRuleCall_1; }
	}
	public class BooleanValueElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.text.ElkGraph.BooleanValue");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final Keyword cTrueKeyword_0 = (Keyword)cAlternatives.eContents().get(0);
		private final Keyword cFalseKeyword_1 = (Keyword)cAlternatives.eContents().get(1);
		
		//BooleanValue ecore::EJavaObject:
		//	'true' | 'false';
		@Override public ParserRule getRule() { return rule; }
		
		//'true' | 'false'
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//'true'
		public Keyword getTrueKeyword_0() { return cTrueKeyword_0; }
		
		//'false'
		public Keyword getFalseKeyword_1() { return cFalseKeyword_1; }
	}
	
	
	private final RootNodeElements pRootNode;
	private final ElkNodeElements pElkNode;
	private final ElkLabelElements pElkLabel;
	private final ElkPortElements pElkPort;
	private final ShapeLayoutElements pShapeLayout;
	private final ElkEdgeElements pElkEdge;
	private final EdgeLayoutElements pEdgeLayout;
	private final ElkSingleEdgeSectionElements pElkSingleEdgeSection;
	private final ElkEdgeSectionElements pElkEdgeSection;
	private final ElkBendPointElements pElkBendPoint;
	private final QualifiedIdElements pQualifiedId;
	private final NumberElements pNumber;
	private final PropertyElements pProperty;
	private final PropertyKeyElements pPropertyKey;
	private final StringValueElements pStringValue;
	private final QualifiedIdValueElements pQualifiedIdValue;
	private final NumberValueElements pNumberValue;
	private final BooleanValueElements pBooleanValue;
	private final TerminalRule tSIGNED_INT;
	private final TerminalRule tFLOAT;
	
	private final Grammar grammar;
	
	private final TerminalsGrammarAccess gaTerminals;

	@Inject
	public ElkGraphGrammarAccess(GrammarProvider grammarProvider,
			TerminalsGrammarAccess gaTerminals) {
		this.grammar = internalFindGrammar(grammarProvider);
		this.gaTerminals = gaTerminals;
		this.pRootNode = new RootNodeElements();
		this.pElkNode = new ElkNodeElements();
		this.pElkLabel = new ElkLabelElements();
		this.pElkPort = new ElkPortElements();
		this.pShapeLayout = new ShapeLayoutElements();
		this.pElkEdge = new ElkEdgeElements();
		this.pEdgeLayout = new EdgeLayoutElements();
		this.pElkSingleEdgeSection = new ElkSingleEdgeSectionElements();
		this.pElkEdgeSection = new ElkEdgeSectionElements();
		this.pElkBendPoint = new ElkBendPointElements();
		this.pQualifiedId = new QualifiedIdElements();
		this.pNumber = new NumberElements();
		this.pProperty = new PropertyElements();
		this.pPropertyKey = new PropertyKeyElements();
		this.pStringValue = new StringValueElements();
		this.pQualifiedIdValue = new QualifiedIdValueElements();
		this.pNumberValue = new NumberValueElements();
		this.pBooleanValue = new BooleanValueElements();
		this.tSIGNED_INT = (TerminalRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.text.ElkGraph.SIGNED_INT");
		this.tFLOAT = (TerminalRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.text.ElkGraph.FLOAT");
	}
	
	protected Grammar internalFindGrammar(GrammarProvider grammarProvider) {
		Grammar grammar = grammarProvider.getGrammar(this);
		while (grammar != null) {
			if ("org.eclipse.elk.graph.text.ElkGraph".equals(grammar.getName())) {
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

	
	//RootNode ElkNode:
	//	{ElkNode} ('graph' identifier=ID)?
	//	ShapeLayout?
	//	properties+=Property* (labels+=ElkLabel | ports+=ElkPort | children+=ElkNode | containedEdges+=ElkEdge)*;
	public RootNodeElements getRootNodeAccess() {
		return pRootNode;
	}
	
	public ParserRule getRootNodeRule() {
		return getRootNodeAccess().getRule();
	}
	
	//ElkNode:
	//	'node' identifier=ID ('{'
	//	ShapeLayout?
	//	properties+=Property* (labels+=ElkLabel | ports+=ElkPort | children+=ElkNode | containedEdges+=ElkEdge)*
	//	'}')?;
	public ElkNodeElements getElkNodeAccess() {
		return pElkNode;
	}
	
	public ParserRule getElkNodeRule() {
		return getElkNodeAccess().getRule();
	}
	
	//ElkLabel:
	//	'label' (identifier=ID ':')? text=STRING ('{'
	//	ShapeLayout?
	//	properties+=Property*
	//	labels+=ElkLabel*
	//	'}')?;
	public ElkLabelElements getElkLabelAccess() {
		return pElkLabel;
	}
	
	public ParserRule getElkLabelRule() {
		return getElkLabelAccess().getRule();
	}
	
	//ElkPort:
	//	'port' identifier=ID ('{'
	//	ShapeLayout?
	//	properties+=Property*
	//	labels+=ElkLabel*
	//	'}')?;
	public ElkPortElements getElkPortAccess() {
		return pElkPort;
	}
	
	public ParserRule getElkPortRule() {
		return getElkPortAccess().getRule();
	}
	
	//fragment ShapeLayout returns ElkShape:
	//	'layout' '[' (('position' ':' x=Number ',' y=Number)?
	//	& ('size' ':' width=Number ',' height=Number)?) ']';
	public ShapeLayoutElements getShapeLayoutAccess() {
		return pShapeLayout;
	}
	
	public ParserRule getShapeLayoutRule() {
		return getShapeLayoutAccess().getRule();
	}
	
	///* SuppressWarnings[BidirectionalReference] */ ElkEdge:
	//	'edge' (identifier=ID ':')?
	//	sources+=[ElkConnectableShape|QualifiedId] (',' sources+=[ElkConnectableShape|QualifiedId])* '->'
	//	targets+=[ElkConnectableShape|QualifiedId] (',' targets+=[ElkConnectableShape|QualifiedId])* ('{'
	//	EdgeLayout?
	//	properties+=Property*
	//	labels+=ElkLabel*
	//	'}')?;
	public ElkEdgeElements getElkEdgeAccess() {
		return pElkEdge;
	}
	
	public ParserRule getElkEdgeRule() {
		return getElkEdgeAccess().getRule();
	}
	
	//fragment EdgeLayout returns ElkEdge:
	//	'layout' '[' (sections+=ElkSingleEdgeSection | sections+=ElkEdgeSection+) ']';
	public EdgeLayoutElements getEdgeLayoutAccess() {
		return pEdgeLayout;
	}
	
	public ParserRule getEdgeLayoutRule() {
		return getEdgeLayoutAccess().getRule();
	}
	
	//ElkSingleEdgeSection ElkEdgeSection:
	//	{ElkEdgeSection} ((('incoming' ':' incomingShape=[ElkConnectableShape|QualifiedId])?
	//	& ('outgoing' ':' outgoingShape=[ElkConnectableShape|QualifiedId])?
	//	& ('start' ':' startX=Number ',' startY=Number)?
	//	& ('end' ':' endX=Number ',' endY=Number)?) ('bends' ':' bendPoints+=ElkBendPoint ('|' bendPoints+=ElkBendPoint)*)?
	//	properties+=Property*);
	public ElkSingleEdgeSectionElements getElkSingleEdgeSectionAccess() {
		return pElkSingleEdgeSection;
	}
	
	public ParserRule getElkSingleEdgeSectionRule() {
		return getElkSingleEdgeSectionAccess().getRule();
	}
	
	///* SuppressWarnings[BidirectionalReference] */ ElkEdgeSection:
	//	'section' identifier=ID ('->' outgoingSections+=[ElkEdgeSection] (',' outgoingSections+=[ElkEdgeSection])*)? '['
	//	((('incoming' ':' incomingShape=[ElkConnectableShape|QualifiedId])?
	//	& ('outgoing' ':' outgoingShape=[ElkConnectableShape|QualifiedId])?
	//	& ('start' ':' startX=Number ',' startY=Number)?
	//	& ('end' ':' endX=Number ',' endY=Number)?) ('bends' ':' bendPoints+=ElkBendPoint ('|' bendPoints+=ElkBendPoint)*)?
	//	properties+=Property*) ']';
	public ElkEdgeSectionElements getElkEdgeSectionAccess() {
		return pElkEdgeSection;
	}
	
	public ParserRule getElkEdgeSectionRule() {
		return getElkEdgeSectionAccess().getRule();
	}
	
	//ElkBendPoint:
	//	x=Number ',' y=Number;
	public ElkBendPointElements getElkBendPointAccess() {
		return pElkBendPoint;
	}
	
	public ParserRule getElkBendPointRule() {
		return getElkBendPointAccess().getRule();
	}
	
	//QualifiedId:
	//	ID ('.' ID)*;
	public QualifiedIdElements getQualifiedIdAccess() {
		return pQualifiedId;
	}
	
	public ParserRule getQualifiedIdRule() {
		return getQualifiedIdAccess().getRule();
	}
	
	//Number ecore::EDouble:
	//	SIGNED_INT | FLOAT;
	public NumberElements getNumberAccess() {
		return pNumber;
	}
	
	public ParserRule getNumberRule() {
		return getNumberAccess().getRule();
	}
	
	//Property ElkPropertyToValueMapEntry:
	//	key=PropertyKey ':' (value=StringValue | value=QualifiedIdValue | value=NumberValue | value=BooleanValue | 'null');
	public PropertyElements getPropertyAccess() {
		return pProperty;
	}
	
	public ParserRule getPropertyRule() {
		return getPropertyAccess().getRule();
	}
	
	//PropertyKey IProperty hidden():
	//	ID ('.' ID)*;
	public PropertyKeyElements getPropertyKeyAccess() {
		return pPropertyKey;
	}
	
	public ParserRule getPropertyKeyRule() {
		return getPropertyKeyAccess().getRule();
	}
	
	//StringValue ecore::EJavaObject:
	//	STRING;
	public StringValueElements getStringValueAccess() {
		return pStringValue;
	}
	
	public ParserRule getStringValueRule() {
		return getStringValueAccess().getRule();
	}
	
	//QualifiedIdValue ecore::EJavaObject:
	//	QualifiedId;
	public QualifiedIdValueElements getQualifiedIdValueAccess() {
		return pQualifiedIdValue;
	}
	
	public ParserRule getQualifiedIdValueRule() {
		return getQualifiedIdValueAccess().getRule();
	}
	
	//NumberValue ecore::EJavaObject:
	//	SIGNED_INT | FLOAT;
	public NumberValueElements getNumberValueAccess() {
		return pNumberValue;
	}
	
	public ParserRule getNumberValueRule() {
		return getNumberValueAccess().getRule();
	}
	
	//BooleanValue ecore::EJavaObject:
	//	'true' | 'false';
	public BooleanValueElements getBooleanValueAccess() {
		return pBooleanValue;
	}
	
	public ParserRule getBooleanValueRule() {
		return getBooleanValueAccess().getRule();
	}
	
	//terminal SIGNED_INT returns ecore::EInt:
	//	('+' | '-')? INT;
	public TerminalRule getSIGNED_INTRule() {
		return tSIGNED_INT;
	}
	
	//terminal FLOAT returns ecore::EDouble:
	//	('+' | '-')? (INT '.' INT | INT ('.' INT)? ('e' | 'E') ('+' | '-')? INT);
	public TerminalRule getFLOATRule() {
		return tFLOAT;
	}
	
	//terminal ID:
	//	'^'? ('a'..'z' | 'A'..'Z' | '_') ('a'..'z' | 'A'..'Z' | '_' | '0'..'9')*;
	public TerminalRule getIDRule() {
		return gaTerminals.getIDRule();
	}
	
	//terminal INT returns ecore::EInt:
	//	'0'..'9'+;
	public TerminalRule getINTRule() {
		return gaTerminals.getINTRule();
	}
	
	//terminal STRING:
	//	'"' ('\\' . | !('\\' | '"'))* '"' |
	//	"'" ('\\' . | !('\\' | "'"))* "'";
	public TerminalRule getSTRINGRule() {
		return gaTerminals.getSTRINGRule();
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
