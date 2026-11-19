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
package org.eclipse.elk.graph.json.text.services;

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
import org.eclipse.xtext.common.services.TerminalsGrammarAccess;
import org.eclipse.xtext.service.AbstractElementFinder.AbstractGrammarElementFinder;
import org.eclipse.xtext.service.GrammarProvider;

@Singleton
public class ElkGraphJsonGrammarAccess extends AbstractGrammarElementFinder {
	
	public class ElkNodeElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.ElkNode");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Action cElkNodeAction_0 = (Action)cGroup.eContents().get(0);
		private final Keyword cLeftCurlyBracketKeyword_1 = (Keyword)cGroup.eContents().get(1);
		private final Group cGroup_2 = (Group)cGroup.eContents().get(2);
		private final RuleCall cNodeElementParserRuleCall_2_0 = (RuleCall)cGroup_2.eContents().get(0);
		private final Group cGroup_2_1 = (Group)cGroup_2.eContents().get(1);
		private final Keyword cCommaKeyword_2_1_0 = (Keyword)cGroup_2_1.eContents().get(0);
		private final RuleCall cNodeElementParserRuleCall_2_1_1 = (RuleCall)cGroup_2_1.eContents().get(1);
		private final Keyword cCommaKeyword_3 = (Keyword)cGroup.eContents().get(3);
		private final Keyword cRightCurlyBracketKeyword_4 = (Keyword)cGroup.eContents().get(4);
		
		///* - - Nodes - - */ ElkNode:
		//	{ElkNode} '{' (NodeElement (',' NodeElement)*)? ','? '}';
		@Override public ParserRule getRule() { return rule; }
		
		//{ElkNode} '{' (NodeElement (',' NodeElement)*)? ','? '}'
		public Group getGroup() { return cGroup; }
		
		//{ElkNode}
		public Action getElkNodeAction_0() { return cElkNodeAction_0; }
		
		//'{'
		public Keyword getLeftCurlyBracketKeyword_1() { return cLeftCurlyBracketKeyword_1; }
		
		//(NodeElement (',' NodeElement)*)?
		public Group getGroup_2() { return cGroup_2; }
		
		//NodeElement
		public RuleCall getNodeElementParserRuleCall_2_0() { return cNodeElementParserRuleCall_2_0; }
		
		//(',' NodeElement)*
		public Group getGroup_2_1() { return cGroup_2_1; }
		
		//','
		public Keyword getCommaKeyword_2_1_0() { return cCommaKeyword_2_1_0; }
		
		//NodeElement
		public RuleCall getNodeElementParserRuleCall_2_1_1() { return cNodeElementParserRuleCall_2_1_1; }
		
		//','?
		public Keyword getCommaKeyword_3() { return cCommaKeyword_3; }
		
		//'}'
		public Keyword getRightCurlyBracketKeyword_4() { return cRightCurlyBracketKeyword_4; }
	}
	public class NodeElementElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.NodeElement");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final RuleCall cElkIdParserRuleCall_0 = (RuleCall)cAlternatives.eContents().get(0);
		private final RuleCall cShapeElementParserRuleCall_1 = (RuleCall)cAlternatives.eContents().get(1);
		private final Group cGroup_2 = (Group)cAlternatives.eContents().get(2);
		private final RuleCall cKeyChildrenParserRuleCall_2_0 = (RuleCall)cGroup_2.eContents().get(0);
		private final Keyword cColonKeyword_2_1 = (Keyword)cGroup_2.eContents().get(1);
		private final RuleCall cElkNodeChildrenParserRuleCall_2_2 = (RuleCall)cGroup_2.eContents().get(2);
		private final Group cGroup_3 = (Group)cAlternatives.eContents().get(3);
		private final RuleCall cKeyPortsParserRuleCall_3_0 = (RuleCall)cGroup_3.eContents().get(0);
		private final Keyword cColonKeyword_3_1 = (Keyword)cGroup_3.eContents().get(1);
		private final RuleCall cElkNodePortsParserRuleCall_3_2 = (RuleCall)cGroup_3.eContents().get(2);
		private final Group cGroup_4 = (Group)cAlternatives.eContents().get(4);
		private final RuleCall cKeyLabelsParserRuleCall_4_0 = (RuleCall)cGroup_4.eContents().get(0);
		private final Keyword cColonKeyword_4_1 = (Keyword)cGroup_4.eContents().get(1);
		private final RuleCall cElkGraphElementLabelsParserRuleCall_4_2 = (RuleCall)cGroup_4.eContents().get(2);
		private final Group cGroup_5 = (Group)cAlternatives.eContents().get(5);
		private final RuleCall cKeyEdgesParserRuleCall_5_0 = (RuleCall)cGroup_5.eContents().get(0);
		private final Keyword cColonKeyword_5_1 = (Keyword)cGroup_5.eContents().get(1);
		private final RuleCall cElkNodeEdgesParserRuleCall_5_2 = (RuleCall)cGroup_5.eContents().get(2);
		private final Group cGroup_6 = (Group)cAlternatives.eContents().get(6);
		private final RuleCall cKeyLayoutOptionsParserRuleCall_6_0 = (RuleCall)cGroup_6.eContents().get(0);
		private final Keyword cColonKeyword_6_1 = (Keyword)cGroup_6.eContents().get(1);
		private final RuleCall cElkGraphElementPropertiesParserRuleCall_6_2 = (RuleCall)cGroup_6.eContents().get(2);
		private final RuleCall cJsonMemberParserRuleCall_7 = (RuleCall)cAlternatives.eContents().get(7);
		
		//fragment NodeElement returns ElkNode:
		//	ElkId
		//	| ShapeElement
		//	| KeyChildren ':' ElkNodeChildren
		//	| KeyPorts ':' ElkNodePorts
		//	| KeyLabels ':' ElkGraphElementLabels
		//	| KeyEdges ':' ElkNodeEdges
		//	| KeyLayoutOptions ':' ElkGraphElementProperties
		//	| JsonMember;
		@Override public ParserRule getRule() { return rule; }
		
		//ElkId | ShapeElement | KeyChildren ':' ElkNodeChildren | KeyPorts ':' ElkNodePorts | KeyLabels ':' ElkGraphElementLabels
		//| KeyEdges ':' ElkNodeEdges | KeyLayoutOptions ':' ElkGraphElementProperties | JsonMember
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//ElkId
		public RuleCall getElkIdParserRuleCall_0() { return cElkIdParserRuleCall_0; }
		
		//ShapeElement
		public RuleCall getShapeElementParserRuleCall_1() { return cShapeElementParserRuleCall_1; }
		
		//KeyChildren ':' ElkNodeChildren
		public Group getGroup_2() { return cGroup_2; }
		
		//KeyChildren
		public RuleCall getKeyChildrenParserRuleCall_2_0() { return cKeyChildrenParserRuleCall_2_0; }
		
		//':'
		public Keyword getColonKeyword_2_1() { return cColonKeyword_2_1; }
		
		//ElkNodeChildren
		public RuleCall getElkNodeChildrenParserRuleCall_2_2() { return cElkNodeChildrenParserRuleCall_2_2; }
		
		//KeyPorts ':' ElkNodePorts
		public Group getGroup_3() { return cGroup_3; }
		
		//KeyPorts
		public RuleCall getKeyPortsParserRuleCall_3_0() { return cKeyPortsParserRuleCall_3_0; }
		
		//':'
		public Keyword getColonKeyword_3_1() { return cColonKeyword_3_1; }
		
		//ElkNodePorts
		public RuleCall getElkNodePortsParserRuleCall_3_2() { return cElkNodePortsParserRuleCall_3_2; }
		
		//KeyLabels ':' ElkGraphElementLabels
		public Group getGroup_4() { return cGroup_4; }
		
		//KeyLabels
		public RuleCall getKeyLabelsParserRuleCall_4_0() { return cKeyLabelsParserRuleCall_4_0; }
		
		//':'
		public Keyword getColonKeyword_4_1() { return cColonKeyword_4_1; }
		
		//ElkGraphElementLabels
		public RuleCall getElkGraphElementLabelsParserRuleCall_4_2() { return cElkGraphElementLabelsParserRuleCall_4_2; }
		
		//KeyEdges ':' ElkNodeEdges
		public Group getGroup_5() { return cGroup_5; }
		
		//KeyEdges
		public RuleCall getKeyEdgesParserRuleCall_5_0() { return cKeyEdgesParserRuleCall_5_0; }
		
		//':'
		public Keyword getColonKeyword_5_1() { return cColonKeyword_5_1; }
		
		//ElkNodeEdges
		public RuleCall getElkNodeEdgesParserRuleCall_5_2() { return cElkNodeEdgesParserRuleCall_5_2; }
		
		//KeyLayoutOptions ':' ElkGraphElementProperties
		public Group getGroup_6() { return cGroup_6; }
		
		//KeyLayoutOptions
		public RuleCall getKeyLayoutOptionsParserRuleCall_6_0() { return cKeyLayoutOptionsParserRuleCall_6_0; }
		
		//':'
		public Keyword getColonKeyword_6_1() { return cColonKeyword_6_1; }
		
		//ElkGraphElementProperties
		public RuleCall getElkGraphElementPropertiesParserRuleCall_6_2() { return cElkGraphElementPropertiesParserRuleCall_6_2; }
		
		//JsonMember
		public RuleCall getJsonMemberParserRuleCall_7() { return cJsonMemberParserRuleCall_7; }
	}
	public class ElkPortElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.ElkPort");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Keyword cLeftCurlyBracketKeyword_0 = (Keyword)cGroup.eContents().get(0);
		private final RuleCall cPortElementParserRuleCall_1 = (RuleCall)cGroup.eContents().get(1);
		private final Group cGroup_2 = (Group)cGroup.eContents().get(2);
		private final Keyword cCommaKeyword_2_0 = (Keyword)cGroup_2.eContents().get(0);
		private final RuleCall cPortElementParserRuleCall_2_1 = (RuleCall)cGroup_2.eContents().get(1);
		private final Keyword cCommaKeyword_3 = (Keyword)cGroup.eContents().get(3);
		private final Keyword cRightCurlyBracketKeyword_4 = (Keyword)cGroup.eContents().get(4);
		
		///* - - Ports - - */ ElkPort:
		//	'{' PortElement (',' PortElement)* ','? '}';
		@Override public ParserRule getRule() { return rule; }
		
		//'{' PortElement (',' PortElement)* ','? '}'
		public Group getGroup() { return cGroup; }
		
		//'{'
		public Keyword getLeftCurlyBracketKeyword_0() { return cLeftCurlyBracketKeyword_0; }
		
		//PortElement
		public RuleCall getPortElementParserRuleCall_1() { return cPortElementParserRuleCall_1; }
		
		//(',' PortElement)*
		public Group getGroup_2() { return cGroup_2; }
		
		//','
		public Keyword getCommaKeyword_2_0() { return cCommaKeyword_2_0; }
		
		//PortElement
		public RuleCall getPortElementParserRuleCall_2_1() { return cPortElementParserRuleCall_2_1; }
		
		//','?
		public Keyword getCommaKeyword_3() { return cCommaKeyword_3; }
		
		//'}'
		public Keyword getRightCurlyBracketKeyword_4() { return cRightCurlyBracketKeyword_4; }
	}
	public class PortElementElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.PortElement");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final RuleCall cElkIdParserRuleCall_0 = (RuleCall)cAlternatives.eContents().get(0);
		private final RuleCall cShapeElementParserRuleCall_1 = (RuleCall)cAlternatives.eContents().get(1);
		private final Group cGroup_2 = (Group)cAlternatives.eContents().get(2);
		private final RuleCall cKeyLabelsParserRuleCall_2_0 = (RuleCall)cGroup_2.eContents().get(0);
		private final Keyword cColonKeyword_2_1 = (Keyword)cGroup_2.eContents().get(1);
		private final RuleCall cElkGraphElementLabelsParserRuleCall_2_2 = (RuleCall)cGroup_2.eContents().get(2);
		private final Group cGroup_3 = (Group)cAlternatives.eContents().get(3);
		private final RuleCall cKeyLayoutOptionsParserRuleCall_3_0 = (RuleCall)cGroup_3.eContents().get(0);
		private final Keyword cColonKeyword_3_1 = (Keyword)cGroup_3.eContents().get(1);
		private final RuleCall cElkGraphElementPropertiesParserRuleCall_3_2 = (RuleCall)cGroup_3.eContents().get(2);
		private final RuleCall cJsonMemberParserRuleCall_4 = (RuleCall)cAlternatives.eContents().get(4);
		
		//fragment PortElement returns ElkPort:
		//	ElkId
		//	| ShapeElement
		//	| KeyLabels ':' ElkGraphElementLabels
		//	| KeyLayoutOptions ':' ElkGraphElementProperties
		//	| JsonMember;
		@Override public ParserRule getRule() { return rule; }
		
		//ElkId | ShapeElement | KeyLabels ':' ElkGraphElementLabels | KeyLayoutOptions ':' ElkGraphElementProperties | JsonMember
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//ElkId
		public RuleCall getElkIdParserRuleCall_0() { return cElkIdParserRuleCall_0; }
		
		//ShapeElement
		public RuleCall getShapeElementParserRuleCall_1() { return cShapeElementParserRuleCall_1; }
		
		//KeyLabels ':' ElkGraphElementLabels
		public Group getGroup_2() { return cGroup_2; }
		
		//KeyLabels
		public RuleCall getKeyLabelsParserRuleCall_2_0() { return cKeyLabelsParserRuleCall_2_0; }
		
		//':'
		public Keyword getColonKeyword_2_1() { return cColonKeyword_2_1; }
		
		//ElkGraphElementLabels
		public RuleCall getElkGraphElementLabelsParserRuleCall_2_2() { return cElkGraphElementLabelsParserRuleCall_2_2; }
		
		//KeyLayoutOptions ':' ElkGraphElementProperties
		public Group getGroup_3() { return cGroup_3; }
		
		//KeyLayoutOptions
		public RuleCall getKeyLayoutOptionsParserRuleCall_3_0() { return cKeyLayoutOptionsParserRuleCall_3_0; }
		
		//':'
		public Keyword getColonKeyword_3_1() { return cColonKeyword_3_1; }
		
		//ElkGraphElementProperties
		public RuleCall getElkGraphElementPropertiesParserRuleCall_3_2() { return cElkGraphElementPropertiesParserRuleCall_3_2; }
		
		//JsonMember
		public RuleCall getJsonMemberParserRuleCall_4() { return cJsonMemberParserRuleCall_4; }
	}
	public class ElkLabelElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.ElkLabel");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Keyword cLeftCurlyBracketKeyword_0 = (Keyword)cGroup.eContents().get(0);
		private final RuleCall cLabelElementParserRuleCall_1 = (RuleCall)cGroup.eContents().get(1);
		private final Group cGroup_2 = (Group)cGroup.eContents().get(2);
		private final Keyword cCommaKeyword_2_0 = (Keyword)cGroup_2.eContents().get(0);
		private final RuleCall cLabelElementParserRuleCall_2_1 = (RuleCall)cGroup_2.eContents().get(1);
		private final Keyword cCommaKeyword_3 = (Keyword)cGroup.eContents().get(3);
		private final Keyword cRightCurlyBracketKeyword_4 = (Keyword)cGroup.eContents().get(4);
		
		///* - - Labels - - */ ElkLabel:
		//	'{' LabelElement (',' LabelElement)* ','? '}';
		@Override public ParserRule getRule() { return rule; }
		
		//'{' LabelElement (',' LabelElement)* ','? '}'
		public Group getGroup() { return cGroup; }
		
		//'{'
		public Keyword getLeftCurlyBracketKeyword_0() { return cLeftCurlyBracketKeyword_0; }
		
		//LabelElement
		public RuleCall getLabelElementParserRuleCall_1() { return cLabelElementParserRuleCall_1; }
		
		//(',' LabelElement)*
		public Group getGroup_2() { return cGroup_2; }
		
		//','
		public Keyword getCommaKeyword_2_0() { return cCommaKeyword_2_0; }
		
		//LabelElement
		public RuleCall getLabelElementParserRuleCall_2_1() { return cLabelElementParserRuleCall_2_1; }
		
		//','?
		public Keyword getCommaKeyword_3() { return cCommaKeyword_3; }
		
		//'}'
		public Keyword getRightCurlyBracketKeyword_4() { return cRightCurlyBracketKeyword_4; }
	}
	public class LabelElementElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.LabelElement");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final RuleCall cElkIdParserRuleCall_0 = (RuleCall)cAlternatives.eContents().get(0);
		private final RuleCall cShapeElementParserRuleCall_1 = (RuleCall)cAlternatives.eContents().get(1);
		private final Group cGroup_2 = (Group)cAlternatives.eContents().get(2);
		private final RuleCall cKeyTextParserRuleCall_2_0 = (RuleCall)cGroup_2.eContents().get(0);
		private final Keyword cColonKeyword_2_1 = (Keyword)cGroup_2.eContents().get(1);
		private final Assignment cTextAssignment_2_2 = (Assignment)cGroup_2.eContents().get(2);
		private final RuleCall cTextSTRINGTerminalRuleCall_2_2_0 = (RuleCall)cTextAssignment_2_2.eContents().get(0);
		private final Group cGroup_3 = (Group)cAlternatives.eContents().get(3);
		private final RuleCall cKeyLabelsParserRuleCall_3_0 = (RuleCall)cGroup_3.eContents().get(0);
		private final Keyword cColonKeyword_3_1 = (Keyword)cGroup_3.eContents().get(1);
		private final RuleCall cElkGraphElementLabelsParserRuleCall_3_2 = (RuleCall)cGroup_3.eContents().get(2);
		private final Group cGroup_4 = (Group)cAlternatives.eContents().get(4);
		private final RuleCall cKeyLayoutOptionsParserRuleCall_4_0 = (RuleCall)cGroup_4.eContents().get(0);
		private final Keyword cColonKeyword_4_1 = (Keyword)cGroup_4.eContents().get(1);
		private final RuleCall cElkGraphElementPropertiesParserRuleCall_4_2 = (RuleCall)cGroup_4.eContents().get(2);
		private final RuleCall cJsonMemberParserRuleCall_5 = (RuleCall)cAlternatives.eContents().get(5);
		
		//fragment LabelElement returns ElkLabel:
		//	ElkId
		//	| ShapeElement
		//	| KeyText ':' text=STRING
		//	| KeyLabels ':' ElkGraphElementLabels
		//	| KeyLayoutOptions ':' ElkGraphElementProperties
		//	| JsonMember;
		@Override public ParserRule getRule() { return rule; }
		
		//ElkId | ShapeElement | KeyText ':' text=STRING | KeyLabels ':' ElkGraphElementLabels | KeyLayoutOptions ':'
		//ElkGraphElementProperties | JsonMember
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//ElkId
		public RuleCall getElkIdParserRuleCall_0() { return cElkIdParserRuleCall_0; }
		
		//ShapeElement
		public RuleCall getShapeElementParserRuleCall_1() { return cShapeElementParserRuleCall_1; }
		
		//KeyText ':' text=STRING
		public Group getGroup_2() { return cGroup_2; }
		
		//KeyText
		public RuleCall getKeyTextParserRuleCall_2_0() { return cKeyTextParserRuleCall_2_0; }
		
		//':'
		public Keyword getColonKeyword_2_1() { return cColonKeyword_2_1; }
		
		//text=STRING
		public Assignment getTextAssignment_2_2() { return cTextAssignment_2_2; }
		
		//STRING
		public RuleCall getTextSTRINGTerminalRuleCall_2_2_0() { return cTextSTRINGTerminalRuleCall_2_2_0; }
		
		//KeyLabels ':' ElkGraphElementLabels
		public Group getGroup_3() { return cGroup_3; }
		
		//KeyLabels
		public RuleCall getKeyLabelsParserRuleCall_3_0() { return cKeyLabelsParserRuleCall_3_0; }
		
		//':'
		public Keyword getColonKeyword_3_1() { return cColonKeyword_3_1; }
		
		//ElkGraphElementLabels
		public RuleCall getElkGraphElementLabelsParserRuleCall_3_2() { return cElkGraphElementLabelsParserRuleCall_3_2; }
		
		//KeyLayoutOptions ':' ElkGraphElementProperties
		public Group getGroup_4() { return cGroup_4; }
		
		//KeyLayoutOptions
		public RuleCall getKeyLayoutOptionsParserRuleCall_4_0() { return cKeyLayoutOptionsParserRuleCall_4_0; }
		
		//':'
		public Keyword getColonKeyword_4_1() { return cColonKeyword_4_1; }
		
		//ElkGraphElementProperties
		public RuleCall getElkGraphElementPropertiesParserRuleCall_4_2() { return cElkGraphElementPropertiesParserRuleCall_4_2; }
		
		//JsonMember
		public RuleCall getJsonMemberParserRuleCall_5() { return cJsonMemberParserRuleCall_5; }
	}
	public class ElkEdgeElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.ElkEdge");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Keyword cLeftCurlyBracketKeyword_0 = (Keyword)cGroup.eContents().get(0);
		private final RuleCall cEdgeElementParserRuleCall_1 = (RuleCall)cGroup.eContents().get(1);
		private final Group cGroup_2 = (Group)cGroup.eContents().get(2);
		private final Keyword cCommaKeyword_2_0 = (Keyword)cGroup_2.eContents().get(0);
		private final RuleCall cEdgeElementParserRuleCall_2_1 = (RuleCall)cGroup_2.eContents().get(1);
		private final Keyword cCommaKeyword_3 = (Keyword)cGroup.eContents().get(3);
		private final Keyword cRightCurlyBracketKeyword_4 = (Keyword)cGroup.eContents().get(4);
		
		///* - - Edges - - */ ElkEdge:
		//	'{' EdgeElement (',' EdgeElement)* ','? '}';
		@Override public ParserRule getRule() { return rule; }
		
		//'{' EdgeElement (',' EdgeElement)* ','? '}'
		public Group getGroup() { return cGroup; }
		
		//'{'
		public Keyword getLeftCurlyBracketKeyword_0() { return cLeftCurlyBracketKeyword_0; }
		
		//EdgeElement
		public RuleCall getEdgeElementParserRuleCall_1() { return cEdgeElementParserRuleCall_1; }
		
		//(',' EdgeElement)*
		public Group getGroup_2() { return cGroup_2; }
		
		//','
		public Keyword getCommaKeyword_2_0() { return cCommaKeyword_2_0; }
		
		//EdgeElement
		public RuleCall getEdgeElementParserRuleCall_2_1() { return cEdgeElementParserRuleCall_2_1; }
		
		//','?
		public Keyword getCommaKeyword_3() { return cCommaKeyword_3; }
		
		//'}'
		public Keyword getRightCurlyBracketKeyword_4() { return cRightCurlyBracketKeyword_4; }
	}
	public class EdgeElementElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.EdgeElement");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final RuleCall cElkIdParserRuleCall_0 = (RuleCall)cAlternatives.eContents().get(0);
		private final Group cGroup_1 = (Group)cAlternatives.eContents().get(1);
		private final RuleCall cKeySourcesParserRuleCall_1_0 = (RuleCall)cGroup_1.eContents().get(0);
		private final Keyword cColonKeyword_1_1 = (Keyword)cGroup_1.eContents().get(1);
		private final RuleCall cElkEdgeSourcesParserRuleCall_1_2 = (RuleCall)cGroup_1.eContents().get(2);
		private final Group cGroup_2 = (Group)cAlternatives.eContents().get(2);
		private final RuleCall cKeyTargetsParserRuleCall_2_0 = (RuleCall)cGroup_2.eContents().get(0);
		private final Keyword cColonKeyword_2_1 = (Keyword)cGroup_2.eContents().get(1);
		private final RuleCall cElkEdgeTargetsParserRuleCall_2_2 = (RuleCall)cGroup_2.eContents().get(2);
		private final Group cGroup_3 = (Group)cAlternatives.eContents().get(3);
		private final RuleCall cKeyLabelsParserRuleCall_3_0 = (RuleCall)cGroup_3.eContents().get(0);
		private final Keyword cColonKeyword_3_1 = (Keyword)cGroup_3.eContents().get(1);
		private final RuleCall cElkGraphElementLabelsParserRuleCall_3_2 = (RuleCall)cGroup_3.eContents().get(2);
		private final Group cGroup_4 = (Group)cAlternatives.eContents().get(4);
		private final RuleCall cKeyLayoutOptionsParserRuleCall_4_0 = (RuleCall)cGroup_4.eContents().get(0);
		private final Keyword cColonKeyword_4_1 = (Keyword)cGroup_4.eContents().get(1);
		private final RuleCall cElkGraphElementPropertiesParserRuleCall_4_2 = (RuleCall)cGroup_4.eContents().get(2);
		private final RuleCall cJsonMemberParserRuleCall_5 = (RuleCall)cAlternatives.eContents().get(5);
		
		//fragment EdgeElement returns ElkEdge:
		//	ElkId
		//	// TODO EdgeLayout
		//	| KeySources ':' ElkEdgeSources
		//	| KeyTargets ':' ElkEdgeTargets
		//	| KeyLabels ':' ElkGraphElementLabels
		//	| KeyLayoutOptions ':' ElkGraphElementProperties
		//	| JsonMember;
		@Override public ParserRule getRule() { return rule; }
		
		//ElkId // TODO EdgeLayout
		//| KeySources ':' ElkEdgeSources | KeyTargets ':' ElkEdgeTargets | KeyLabels ':' ElkGraphElementLabels | KeyLayoutOptions
		//':' ElkGraphElementProperties | JsonMember
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//ElkId
		public RuleCall getElkIdParserRuleCall_0() { return cElkIdParserRuleCall_0; }
		
		//KeySources ':' ElkEdgeSources
		public Group getGroup_1() { return cGroup_1; }
		
		//KeySources
		public RuleCall getKeySourcesParserRuleCall_1_0() { return cKeySourcesParserRuleCall_1_0; }
		
		//':'
		public Keyword getColonKeyword_1_1() { return cColonKeyword_1_1; }
		
		//ElkEdgeSources
		public RuleCall getElkEdgeSourcesParserRuleCall_1_2() { return cElkEdgeSourcesParserRuleCall_1_2; }
		
		//KeyTargets ':' ElkEdgeTargets
		public Group getGroup_2() { return cGroup_2; }
		
		//KeyTargets
		public RuleCall getKeyTargetsParserRuleCall_2_0() { return cKeyTargetsParserRuleCall_2_0; }
		
		//':'
		public Keyword getColonKeyword_2_1() { return cColonKeyword_2_1; }
		
		//ElkEdgeTargets
		public RuleCall getElkEdgeTargetsParserRuleCall_2_2() { return cElkEdgeTargetsParserRuleCall_2_2; }
		
		//KeyLabels ':' ElkGraphElementLabels
		public Group getGroup_3() { return cGroup_3; }
		
		//KeyLabels
		public RuleCall getKeyLabelsParserRuleCall_3_0() { return cKeyLabelsParserRuleCall_3_0; }
		
		//':'
		public Keyword getColonKeyword_3_1() { return cColonKeyword_3_1; }
		
		//ElkGraphElementLabels
		public RuleCall getElkGraphElementLabelsParserRuleCall_3_2() { return cElkGraphElementLabelsParserRuleCall_3_2; }
		
		//KeyLayoutOptions ':' ElkGraphElementProperties
		public Group getGroup_4() { return cGroup_4; }
		
		//KeyLayoutOptions
		public RuleCall getKeyLayoutOptionsParserRuleCall_4_0() { return cKeyLayoutOptionsParserRuleCall_4_0; }
		
		//':'
		public Keyword getColonKeyword_4_1() { return cColonKeyword_4_1; }
		
		//ElkGraphElementProperties
		public RuleCall getElkGraphElementPropertiesParserRuleCall_4_2() { return cElkGraphElementPropertiesParserRuleCall_4_2; }
		
		//JsonMember
		public RuleCall getJsonMemberParserRuleCall_5() { return cJsonMemberParserRuleCall_5; }
	}
	public class ElkEdgeSourcesElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.ElkEdgeSources");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Keyword cLeftSquareBracketKeyword_0 = (Keyword)cGroup.eContents().get(0);
		private final Group cGroup_1 = (Group)cGroup.eContents().get(1);
		private final Assignment cSourcesAssignment_1_0 = (Assignment)cGroup_1.eContents().get(0);
		private final CrossReference cSourcesElkConnectableShapeCrossReference_1_0_0 = (CrossReference)cSourcesAssignment_1_0.eContents().get(0);
		private final RuleCall cSourcesElkConnectableShapeSTRINGTerminalRuleCall_1_0_0_1 = (RuleCall)cSourcesElkConnectableShapeCrossReference_1_0_0.eContents().get(1);
		private final Group cGroup_1_1 = (Group)cGroup_1.eContents().get(1);
		private final Keyword cCommaKeyword_1_1_0 = (Keyword)cGroup_1_1.eContents().get(0);
		private final Assignment cSourcesAssignment_1_1_1 = (Assignment)cGroup_1_1.eContents().get(1);
		private final CrossReference cSourcesElkConnectableShapeCrossReference_1_1_1_0 = (CrossReference)cSourcesAssignment_1_1_1.eContents().get(0);
		private final RuleCall cSourcesElkConnectableShapeSTRINGTerminalRuleCall_1_1_1_0_1 = (RuleCall)cSourcesElkConnectableShapeCrossReference_1_1_1_0.eContents().get(1);
		private final Keyword cCommaKeyword_2 = (Keyword)cGroup.eContents().get(2);
		private final Keyword cRightSquareBracketKeyword_3 = (Keyword)cGroup.eContents().get(3);
		
		///* SuppressWarnings[BidirectionalReference] */ fragment ElkEdgeSources returns ElkEdge:
		//	'[' (sources+=[ElkConnectableShape|STRING] (',' sources+=[ElkConnectableShape|STRING])*)? ','? ']';
		@Override public ParserRule getRule() { return rule; }
		
		//'[' (sources+=[ElkConnectableShape|STRING] (',' sources+=[ElkConnectableShape|STRING])*)? ','? ']'
		public Group getGroup() { return cGroup; }
		
		//'['
		public Keyword getLeftSquareBracketKeyword_0() { return cLeftSquareBracketKeyword_0; }
		
		//(sources+=[ElkConnectableShape|STRING] (',' sources+=[ElkConnectableShape|STRING])*)?
		public Group getGroup_1() { return cGroup_1; }
		
		//sources+=[ElkConnectableShape|STRING]
		public Assignment getSourcesAssignment_1_0() { return cSourcesAssignment_1_0; }
		
		//[ElkConnectableShape|STRING]
		public CrossReference getSourcesElkConnectableShapeCrossReference_1_0_0() { return cSourcesElkConnectableShapeCrossReference_1_0_0; }
		
		//STRING
		public RuleCall getSourcesElkConnectableShapeSTRINGTerminalRuleCall_1_0_0_1() { return cSourcesElkConnectableShapeSTRINGTerminalRuleCall_1_0_0_1; }
		
		//(',' sources+=[ElkConnectableShape|STRING])*
		public Group getGroup_1_1() { return cGroup_1_1; }
		
		//','
		public Keyword getCommaKeyword_1_1_0() { return cCommaKeyword_1_1_0; }
		
		//sources+=[ElkConnectableShape|STRING]
		public Assignment getSourcesAssignment_1_1_1() { return cSourcesAssignment_1_1_1; }
		
		//[ElkConnectableShape|STRING]
		public CrossReference getSourcesElkConnectableShapeCrossReference_1_1_1_0() { return cSourcesElkConnectableShapeCrossReference_1_1_1_0; }
		
		//STRING
		public RuleCall getSourcesElkConnectableShapeSTRINGTerminalRuleCall_1_1_1_0_1() { return cSourcesElkConnectableShapeSTRINGTerminalRuleCall_1_1_1_0_1; }
		
		//','?
		public Keyword getCommaKeyword_2() { return cCommaKeyword_2; }
		
		//']'
		public Keyword getRightSquareBracketKeyword_3() { return cRightSquareBracketKeyword_3; }
	}
	public class ElkEdgeTargetsElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.ElkEdgeTargets");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Keyword cLeftSquareBracketKeyword_0 = (Keyword)cGroup.eContents().get(0);
		private final Group cGroup_1 = (Group)cGroup.eContents().get(1);
		private final Assignment cTargetsAssignment_1_0 = (Assignment)cGroup_1.eContents().get(0);
		private final CrossReference cTargetsElkConnectableShapeCrossReference_1_0_0 = (CrossReference)cTargetsAssignment_1_0.eContents().get(0);
		private final RuleCall cTargetsElkConnectableShapeSTRINGTerminalRuleCall_1_0_0_1 = (RuleCall)cTargetsElkConnectableShapeCrossReference_1_0_0.eContents().get(1);
		private final Group cGroup_1_1 = (Group)cGroup_1.eContents().get(1);
		private final Keyword cCommaKeyword_1_1_0 = (Keyword)cGroup_1_1.eContents().get(0);
		private final Assignment cTargetsAssignment_1_1_1 = (Assignment)cGroup_1_1.eContents().get(1);
		private final CrossReference cTargetsElkConnectableShapeCrossReference_1_1_1_0 = (CrossReference)cTargetsAssignment_1_1_1.eContents().get(0);
		private final RuleCall cTargetsElkConnectableShapeSTRINGTerminalRuleCall_1_1_1_0_1 = (RuleCall)cTargetsElkConnectableShapeCrossReference_1_1_1_0.eContents().get(1);
		private final Keyword cCommaKeyword_2 = (Keyword)cGroup.eContents().get(2);
		private final Keyword cRightSquareBracketKeyword_3 = (Keyword)cGroup.eContents().get(3);
		
		///* SuppressWarnings[BidirectionalReference] */ fragment ElkEdgeTargets returns ElkEdge:
		//	'[' (targets+=[ElkConnectableShape|STRING] (',' targets+=[ElkConnectableShape|STRING])*)? ','? ']';
		@Override public ParserRule getRule() { return rule; }
		
		//'[' (targets+=[ElkConnectableShape|STRING] (',' targets+=[ElkConnectableShape|STRING])*)? ','? ']'
		public Group getGroup() { return cGroup; }
		
		//'['
		public Keyword getLeftSquareBracketKeyword_0() { return cLeftSquareBracketKeyword_0; }
		
		//(targets+=[ElkConnectableShape|STRING] (',' targets+=[ElkConnectableShape|STRING])*)?
		public Group getGroup_1() { return cGroup_1; }
		
		//targets+=[ElkConnectableShape|STRING]
		public Assignment getTargetsAssignment_1_0() { return cTargetsAssignment_1_0; }
		
		//[ElkConnectableShape|STRING]
		public CrossReference getTargetsElkConnectableShapeCrossReference_1_0_0() { return cTargetsElkConnectableShapeCrossReference_1_0_0; }
		
		//STRING
		public RuleCall getTargetsElkConnectableShapeSTRINGTerminalRuleCall_1_0_0_1() { return cTargetsElkConnectableShapeSTRINGTerminalRuleCall_1_0_0_1; }
		
		//(',' targets+=[ElkConnectableShape|STRING])*
		public Group getGroup_1_1() { return cGroup_1_1; }
		
		//','
		public Keyword getCommaKeyword_1_1_0() { return cCommaKeyword_1_1_0; }
		
		//targets+=[ElkConnectableShape|STRING]
		public Assignment getTargetsAssignment_1_1_1() { return cTargetsAssignment_1_1_1; }
		
		//[ElkConnectableShape|STRING]
		public CrossReference getTargetsElkConnectableShapeCrossReference_1_1_1_0() { return cTargetsElkConnectableShapeCrossReference_1_1_1_0; }
		
		//STRING
		public RuleCall getTargetsElkConnectableShapeSTRINGTerminalRuleCall_1_1_1_0_1() { return cTargetsElkConnectableShapeSTRINGTerminalRuleCall_1_1_1_0_1; }
		
		//','?
		public Keyword getCommaKeyword_2() { return cCommaKeyword_2; }
		
		//']'
		public Keyword getRightSquareBracketKeyword_3() { return cRightSquareBracketKeyword_3; }
	}
	public class ElkIdElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.ElkId");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final RuleCall cKeyIdParserRuleCall_0 = (RuleCall)cGroup.eContents().get(0);
		private final Keyword cColonKeyword_1 = (Keyword)cGroup.eContents().get(1);
		private final Assignment cIdentifierAssignment_2 = (Assignment)cGroup.eContents().get(2);
		private final RuleCall cIdentifierSTRINGTerminalRuleCall_2_0 = (RuleCall)cIdentifierAssignment_2.eContents().get(0);
		
		///* - - ElkGraphElement fragments - - */ fragment ElkId returns ElkGraphElement:
		//	KeyId ':' identifier=STRING;
		@Override public ParserRule getRule() { return rule; }
		
		//KeyId ':' identifier=STRING
		public Group getGroup() { return cGroup; }
		
		//KeyId
		public RuleCall getKeyIdParserRuleCall_0() { return cKeyIdParserRuleCall_0; }
		
		//':'
		public Keyword getColonKeyword_1() { return cColonKeyword_1; }
		
		//identifier=STRING
		public Assignment getIdentifierAssignment_2() { return cIdentifierAssignment_2; }
		
		//STRING
		public RuleCall getIdentifierSTRINGTerminalRuleCall_2_0() { return cIdentifierSTRINGTerminalRuleCall_2_0; }
	}
	public class ElkNodeChildrenElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.ElkNodeChildren");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Keyword cLeftSquareBracketKeyword_0 = (Keyword)cGroup.eContents().get(0);
		private final Group cGroup_1 = (Group)cGroup.eContents().get(1);
		private final Assignment cChildrenAssignment_1_0 = (Assignment)cGroup_1.eContents().get(0);
		private final RuleCall cChildrenElkNodeParserRuleCall_1_0_0 = (RuleCall)cChildrenAssignment_1_0.eContents().get(0);
		private final Group cGroup_1_1 = (Group)cGroup_1.eContents().get(1);
		private final Keyword cCommaKeyword_1_1_0 = (Keyword)cGroup_1_1.eContents().get(0);
		private final Assignment cChildrenAssignment_1_1_1 = (Assignment)cGroup_1_1.eContents().get(1);
		private final RuleCall cChildrenElkNodeParserRuleCall_1_1_1_0 = (RuleCall)cChildrenAssignment_1_1_1.eContents().get(0);
		private final Keyword cCommaKeyword_2 = (Keyword)cGroup.eContents().get(2);
		private final Keyword cRightSquareBracketKeyword_3 = (Keyword)cGroup.eContents().get(3);
		
		//fragment ElkNodeChildren returns ElkNode:
		//	'[' (children+=ElkNode (',' children+=ElkNode)*)? ','? ']';
		@Override public ParserRule getRule() { return rule; }
		
		//'[' (children+=ElkNode (',' children+=ElkNode)*)? ','? ']'
		public Group getGroup() { return cGroup; }
		
		//'['
		public Keyword getLeftSquareBracketKeyword_0() { return cLeftSquareBracketKeyword_0; }
		
		//(children+=ElkNode (',' children+=ElkNode)*)?
		public Group getGroup_1() { return cGroup_1; }
		
		//children+=ElkNode
		public Assignment getChildrenAssignment_1_0() { return cChildrenAssignment_1_0; }
		
		//ElkNode
		public RuleCall getChildrenElkNodeParserRuleCall_1_0_0() { return cChildrenElkNodeParserRuleCall_1_0_0; }
		
		//(',' children+=ElkNode)*
		public Group getGroup_1_1() { return cGroup_1_1; }
		
		//','
		public Keyword getCommaKeyword_1_1_0() { return cCommaKeyword_1_1_0; }
		
		//children+=ElkNode
		public Assignment getChildrenAssignment_1_1_1() { return cChildrenAssignment_1_1_1; }
		
		//ElkNode
		public RuleCall getChildrenElkNodeParserRuleCall_1_1_1_0() { return cChildrenElkNodeParserRuleCall_1_1_1_0; }
		
		//','?
		public Keyword getCommaKeyword_2() { return cCommaKeyword_2; }
		
		//']'
		public Keyword getRightSquareBracketKeyword_3() { return cRightSquareBracketKeyword_3; }
	}
	public class ElkNodePortsElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.ElkNodePorts");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Keyword cLeftSquareBracketKeyword_0 = (Keyword)cGroup.eContents().get(0);
		private final Group cGroup_1 = (Group)cGroup.eContents().get(1);
		private final Assignment cPortsAssignment_1_0 = (Assignment)cGroup_1.eContents().get(0);
		private final RuleCall cPortsElkPortParserRuleCall_1_0_0 = (RuleCall)cPortsAssignment_1_0.eContents().get(0);
		private final Group cGroup_1_1 = (Group)cGroup_1.eContents().get(1);
		private final Keyword cCommaKeyword_1_1_0 = (Keyword)cGroup_1_1.eContents().get(0);
		private final Assignment cPortsAssignment_1_1_1 = (Assignment)cGroup_1_1.eContents().get(1);
		private final RuleCall cPortsElkPortParserRuleCall_1_1_1_0 = (RuleCall)cPortsAssignment_1_1_1.eContents().get(0);
		private final Keyword cCommaKeyword_2 = (Keyword)cGroup.eContents().get(2);
		private final Keyword cRightSquareBracketKeyword_3 = (Keyword)cGroup.eContents().get(3);
		
		//fragment ElkNodePorts returns ElkNode:
		//	'[' (ports+=ElkPort (',' ports+=ElkPort)*)? ','? ']';
		@Override public ParserRule getRule() { return rule; }
		
		//'[' (ports+=ElkPort (',' ports+=ElkPort)*)? ','? ']'
		public Group getGroup() { return cGroup; }
		
		//'['
		public Keyword getLeftSquareBracketKeyword_0() { return cLeftSquareBracketKeyword_0; }
		
		//(ports+=ElkPort (',' ports+=ElkPort)*)?
		public Group getGroup_1() { return cGroup_1; }
		
		//ports+=ElkPort
		public Assignment getPortsAssignment_1_0() { return cPortsAssignment_1_0; }
		
		//ElkPort
		public RuleCall getPortsElkPortParserRuleCall_1_0_0() { return cPortsElkPortParserRuleCall_1_0_0; }
		
		//(',' ports+=ElkPort)*
		public Group getGroup_1_1() { return cGroup_1_1; }
		
		//','
		public Keyword getCommaKeyword_1_1_0() { return cCommaKeyword_1_1_0; }
		
		//ports+=ElkPort
		public Assignment getPortsAssignment_1_1_1() { return cPortsAssignment_1_1_1; }
		
		//ElkPort
		public RuleCall getPortsElkPortParserRuleCall_1_1_1_0() { return cPortsElkPortParserRuleCall_1_1_1_0; }
		
		//','?
		public Keyword getCommaKeyword_2() { return cCommaKeyword_2; }
		
		//']'
		public Keyword getRightSquareBracketKeyword_3() { return cRightSquareBracketKeyword_3; }
	}
	public class ElkNodeEdgesElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.ElkNodeEdges");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Keyword cLeftSquareBracketKeyword_0 = (Keyword)cGroup.eContents().get(0);
		private final Group cGroup_1 = (Group)cGroup.eContents().get(1);
		private final Assignment cContainedEdgesAssignment_1_0 = (Assignment)cGroup_1.eContents().get(0);
		private final RuleCall cContainedEdgesElkEdgeParserRuleCall_1_0_0 = (RuleCall)cContainedEdgesAssignment_1_0.eContents().get(0);
		private final Group cGroup_1_1 = (Group)cGroup_1.eContents().get(1);
		private final Keyword cCommaKeyword_1_1_0 = (Keyword)cGroup_1_1.eContents().get(0);
		private final Assignment cContainedEdgesAssignment_1_1_1 = (Assignment)cGroup_1_1.eContents().get(1);
		private final RuleCall cContainedEdgesElkEdgeParserRuleCall_1_1_1_0 = (RuleCall)cContainedEdgesAssignment_1_1_1.eContents().get(0);
		private final Keyword cCommaKeyword_2 = (Keyword)cGroup.eContents().get(2);
		private final Keyword cRightSquareBracketKeyword_3 = (Keyword)cGroup.eContents().get(3);
		
		//fragment ElkNodeEdges returns ElkNode:
		//	'[' (containedEdges+=ElkEdge (',' containedEdges+=ElkEdge)*)? ','? ']';
		@Override public ParserRule getRule() { return rule; }
		
		//'[' (containedEdges+=ElkEdge (',' containedEdges+=ElkEdge)*)? ','? ']'
		public Group getGroup() { return cGroup; }
		
		//'['
		public Keyword getLeftSquareBracketKeyword_0() { return cLeftSquareBracketKeyword_0; }
		
		//(containedEdges+=ElkEdge (',' containedEdges+=ElkEdge)*)?
		public Group getGroup_1() { return cGroup_1; }
		
		//containedEdges+=ElkEdge
		public Assignment getContainedEdgesAssignment_1_0() { return cContainedEdgesAssignment_1_0; }
		
		//ElkEdge
		public RuleCall getContainedEdgesElkEdgeParserRuleCall_1_0_0() { return cContainedEdgesElkEdgeParserRuleCall_1_0_0; }
		
		//(',' containedEdges+=ElkEdge)*
		public Group getGroup_1_1() { return cGroup_1_1; }
		
		//','
		public Keyword getCommaKeyword_1_1_0() { return cCommaKeyword_1_1_0; }
		
		//containedEdges+=ElkEdge
		public Assignment getContainedEdgesAssignment_1_1_1() { return cContainedEdgesAssignment_1_1_1; }
		
		//ElkEdge
		public RuleCall getContainedEdgesElkEdgeParserRuleCall_1_1_1_0() { return cContainedEdgesElkEdgeParserRuleCall_1_1_1_0; }
		
		//','?
		public Keyword getCommaKeyword_2() { return cCommaKeyword_2; }
		
		//']'
		public Keyword getRightSquareBracketKeyword_3() { return cRightSquareBracketKeyword_3; }
	}
	public class ElkGraphElementLabelsElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.ElkGraphElementLabels");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Keyword cLeftSquareBracketKeyword_0 = (Keyword)cGroup.eContents().get(0);
		private final Group cGroup_1 = (Group)cGroup.eContents().get(1);
		private final Assignment cLabelsAssignment_1_0 = (Assignment)cGroup_1.eContents().get(0);
		private final RuleCall cLabelsElkLabelParserRuleCall_1_0_0 = (RuleCall)cLabelsAssignment_1_0.eContents().get(0);
		private final Group cGroup_1_1 = (Group)cGroup_1.eContents().get(1);
		private final Keyword cCommaKeyword_1_1_0 = (Keyword)cGroup_1_1.eContents().get(0);
		private final Assignment cLabelsAssignment_1_1_1 = (Assignment)cGroup_1_1.eContents().get(1);
		private final RuleCall cLabelsElkLabelParserRuleCall_1_1_1_0 = (RuleCall)cLabelsAssignment_1_1_1.eContents().get(0);
		private final Keyword cCommaKeyword_2 = (Keyword)cGroup.eContents().get(2);
		private final Keyword cRightSquareBracketKeyword_3 = (Keyword)cGroup.eContents().get(3);
		
		//fragment ElkGraphElementLabels returns ElkGraphElement:
		//	'[' (labels+=ElkLabel (',' labels+=ElkLabel)*)? ','? ']';
		@Override public ParserRule getRule() { return rule; }
		
		//'[' (labels+=ElkLabel (',' labels+=ElkLabel)*)? ','? ']'
		public Group getGroup() { return cGroup; }
		
		//'['
		public Keyword getLeftSquareBracketKeyword_0() { return cLeftSquareBracketKeyword_0; }
		
		//(labels+=ElkLabel (',' labels+=ElkLabel)*)?
		public Group getGroup_1() { return cGroup_1; }
		
		//labels+=ElkLabel
		public Assignment getLabelsAssignment_1_0() { return cLabelsAssignment_1_0; }
		
		//ElkLabel
		public RuleCall getLabelsElkLabelParserRuleCall_1_0_0() { return cLabelsElkLabelParserRuleCall_1_0_0; }
		
		//(',' labels+=ElkLabel)*
		public Group getGroup_1_1() { return cGroup_1_1; }
		
		//','
		public Keyword getCommaKeyword_1_1_0() { return cCommaKeyword_1_1_0; }
		
		//labels+=ElkLabel
		public Assignment getLabelsAssignment_1_1_1() { return cLabelsAssignment_1_1_1; }
		
		//ElkLabel
		public RuleCall getLabelsElkLabelParserRuleCall_1_1_1_0() { return cLabelsElkLabelParserRuleCall_1_1_1_0; }
		
		//','?
		public Keyword getCommaKeyword_2() { return cCommaKeyword_2; }
		
		//']'
		public Keyword getRightSquareBracketKeyword_3() { return cRightSquareBracketKeyword_3; }
	}
	public class ElkGraphElementPropertiesElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.ElkGraphElementProperties");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Keyword cLeftCurlyBracketKeyword_0 = (Keyword)cGroup.eContents().get(0);
		private final Group cGroup_1 = (Group)cGroup.eContents().get(1);
		private final Assignment cPropertiesAssignment_1_0 = (Assignment)cGroup_1.eContents().get(0);
		private final RuleCall cPropertiesPropertyParserRuleCall_1_0_0 = (RuleCall)cPropertiesAssignment_1_0.eContents().get(0);
		private final Group cGroup_1_1 = (Group)cGroup_1.eContents().get(1);
		private final Keyword cCommaKeyword_1_1_0 = (Keyword)cGroup_1_1.eContents().get(0);
		private final Assignment cPropertiesAssignment_1_1_1 = (Assignment)cGroup_1_1.eContents().get(1);
		private final RuleCall cPropertiesPropertyParserRuleCall_1_1_1_0 = (RuleCall)cPropertiesAssignment_1_1_1.eContents().get(0);
		private final Keyword cCommaKeyword_2 = (Keyword)cGroup.eContents().get(2);
		private final Keyword cRightCurlyBracketKeyword_3 = (Keyword)cGroup.eContents().get(3);
		
		//fragment ElkGraphElementProperties returns ElkGraphElement:
		//	'{' (properties+=Property (',' properties+=Property)*)? ','? '}';
		@Override public ParserRule getRule() { return rule; }
		
		//'{' (properties+=Property (',' properties+=Property)*)? ','? '}'
		public Group getGroup() { return cGroup; }
		
		//'{'
		public Keyword getLeftCurlyBracketKeyword_0() { return cLeftCurlyBracketKeyword_0; }
		
		//(properties+=Property (',' properties+=Property)*)?
		public Group getGroup_1() { return cGroup_1; }
		
		//properties+=Property
		public Assignment getPropertiesAssignment_1_0() { return cPropertiesAssignment_1_0; }
		
		//Property
		public RuleCall getPropertiesPropertyParserRuleCall_1_0_0() { return cPropertiesPropertyParserRuleCall_1_0_0; }
		
		//(',' properties+=Property)*
		public Group getGroup_1_1() { return cGroup_1_1; }
		
		//','
		public Keyword getCommaKeyword_1_1_0() { return cCommaKeyword_1_1_0; }
		
		//properties+=Property
		public Assignment getPropertiesAssignment_1_1_1() { return cPropertiesAssignment_1_1_1; }
		
		//Property
		public RuleCall getPropertiesPropertyParserRuleCall_1_1_1_0() { return cPropertiesPropertyParserRuleCall_1_1_1_0; }
		
		//','?
		public Keyword getCommaKeyword_2() { return cCommaKeyword_2; }
		
		//'}'
		public Keyword getRightCurlyBracketKeyword_3() { return cRightCurlyBracketKeyword_3; }
	}
	public class ShapeElementElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.ShapeElement");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final Group cGroup_0 = (Group)cAlternatives.eContents().get(0);
		private final RuleCall cKeyXParserRuleCall_0_0 = (RuleCall)cGroup_0.eContents().get(0);
		private final Keyword cColonKeyword_0_1 = (Keyword)cGroup_0.eContents().get(1);
		private final Assignment cXAssignment_0_2 = (Assignment)cGroup_0.eContents().get(2);
		private final RuleCall cXNumberParserRuleCall_0_2_0 = (RuleCall)cXAssignment_0_2.eContents().get(0);
		private final Group cGroup_1 = (Group)cAlternatives.eContents().get(1);
		private final RuleCall cKeyYParserRuleCall_1_0 = (RuleCall)cGroup_1.eContents().get(0);
		private final Keyword cColonKeyword_1_1 = (Keyword)cGroup_1.eContents().get(1);
		private final Assignment cYAssignment_1_2 = (Assignment)cGroup_1.eContents().get(2);
		private final RuleCall cYNumberParserRuleCall_1_2_0 = (RuleCall)cYAssignment_1_2.eContents().get(0);
		private final Group cGroup_2 = (Group)cAlternatives.eContents().get(2);
		private final RuleCall cKeyWidthParserRuleCall_2_0 = (RuleCall)cGroup_2.eContents().get(0);
		private final Keyword cColonKeyword_2_1 = (Keyword)cGroup_2.eContents().get(1);
		private final Assignment cWidthAssignment_2_2 = (Assignment)cGroup_2.eContents().get(2);
		private final RuleCall cWidthNumberParserRuleCall_2_2_0 = (RuleCall)cWidthAssignment_2_2.eContents().get(0);
		private final Group cGroup_3 = (Group)cAlternatives.eContents().get(3);
		private final RuleCall cKeyHeightParserRuleCall_3_0 = (RuleCall)cGroup_3.eContents().get(0);
		private final Keyword cColonKeyword_3_1 = (Keyword)cGroup_3.eContents().get(1);
		private final Assignment cHeightAssignment_3_2 = (Assignment)cGroup_3.eContents().get(2);
		private final RuleCall cHeightNumberParserRuleCall_3_2_0 = (RuleCall)cHeightAssignment_3_2.eContents().get(0);
		
		///* - - Shape - - */ fragment ShapeElement returns ElkShape:
		//	KeyX ':' x=Number
		//	| KeyY ':' y=Number
		//	| KeyWidth ':' width=Number
		//	| KeyHeight ':' height=Number;
		@Override public ParserRule getRule() { return rule; }
		
		//KeyX ':' x=Number | KeyY ':' y=Number | KeyWidth ':' width=Number | KeyHeight ':' height=Number
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//KeyX ':' x=Number
		public Group getGroup_0() { return cGroup_0; }
		
		//KeyX
		public RuleCall getKeyXParserRuleCall_0_0() { return cKeyXParserRuleCall_0_0; }
		
		//':'
		public Keyword getColonKeyword_0_1() { return cColonKeyword_0_1; }
		
		//x=Number
		public Assignment getXAssignment_0_2() { return cXAssignment_0_2; }
		
		//Number
		public RuleCall getXNumberParserRuleCall_0_2_0() { return cXNumberParserRuleCall_0_2_0; }
		
		//KeyY ':' y=Number
		public Group getGroup_1() { return cGroup_1; }
		
		//KeyY
		public RuleCall getKeyYParserRuleCall_1_0() { return cKeyYParserRuleCall_1_0; }
		
		//':'
		public Keyword getColonKeyword_1_1() { return cColonKeyword_1_1; }
		
		//y=Number
		public Assignment getYAssignment_1_2() { return cYAssignment_1_2; }
		
		//Number
		public RuleCall getYNumberParserRuleCall_1_2_0() { return cYNumberParserRuleCall_1_2_0; }
		
		//KeyWidth ':' width=Number
		public Group getGroup_2() { return cGroup_2; }
		
		//KeyWidth
		public RuleCall getKeyWidthParserRuleCall_2_0() { return cKeyWidthParserRuleCall_2_0; }
		
		//':'
		public Keyword getColonKeyword_2_1() { return cColonKeyword_2_1; }
		
		//width=Number
		public Assignment getWidthAssignment_2_2() { return cWidthAssignment_2_2; }
		
		//Number
		public RuleCall getWidthNumberParserRuleCall_2_2_0() { return cWidthNumberParserRuleCall_2_2_0; }
		
		//KeyHeight ':' height=Number
		public Group getGroup_3() { return cGroup_3; }
		
		//KeyHeight
		public RuleCall getKeyHeightParserRuleCall_3_0() { return cKeyHeightParserRuleCall_3_0; }
		
		//':'
		public Keyword getColonKeyword_3_1() { return cColonKeyword_3_1; }
		
		//height=Number
		public Assignment getHeightAssignment_3_2() { return cHeightAssignment_3_2; }
		
		//Number
		public RuleCall getHeightNumberParserRuleCall_3_2_0() { return cHeightNumberParserRuleCall_3_2_0; }
	}
	public class PropertyElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.Property");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Assignment cKeyAssignment_0 = (Assignment)cGroup.eContents().get(0);
		private final RuleCall cKeyPropertyKeyParserRuleCall_0_0 = (RuleCall)cKeyAssignment_0.eContents().get(0);
		private final Keyword cColonKeyword_1 = (Keyword)cGroup.eContents().get(1);
		private final Alternatives cAlternatives_2 = (Alternatives)cGroup.eContents().get(2);
		private final Assignment cValueAssignment_2_0 = (Assignment)cAlternatives_2.eContents().get(0);
		private final RuleCall cValueStringValueParserRuleCall_2_0_0 = (RuleCall)cValueAssignment_2_0.eContents().get(0);
		private final Assignment cValueAssignment_2_1 = (Assignment)cAlternatives_2.eContents().get(1);
		private final RuleCall cValueNumberValueParserRuleCall_2_1_0 = (RuleCall)cValueAssignment_2_1.eContents().get(0);
		private final Assignment cValueAssignment_2_2 = (Assignment)cAlternatives_2.eContents().get(2);
		private final RuleCall cValueBooleanValueParserRuleCall_2_2_0 = (RuleCall)cValueAssignment_2_2.eContents().get(0);
		
		//Property ElkPropertyToValueMapEntry:
		//	key=PropertyKey ':' (value=StringValue | value=NumberValue | value=BooleanValue);
		@Override public ParserRule getRule() { return rule; }
		
		//key=PropertyKey ':' (value=StringValue | value=NumberValue | value=BooleanValue)
		public Group getGroup() { return cGroup; }
		
		//key=PropertyKey
		public Assignment getKeyAssignment_0() { return cKeyAssignment_0; }
		
		//PropertyKey
		public RuleCall getKeyPropertyKeyParserRuleCall_0_0() { return cKeyPropertyKeyParserRuleCall_0_0; }
		
		//':'
		public Keyword getColonKeyword_1() { return cColonKeyword_1; }
		
		//(value=StringValue | value=NumberValue | value=BooleanValue)
		public Alternatives getAlternatives_2() { return cAlternatives_2; }
		
		//value=StringValue
		public Assignment getValueAssignment_2_0() { return cValueAssignment_2_0; }
		
		//StringValue
		public RuleCall getValueStringValueParserRuleCall_2_0_0() { return cValueStringValueParserRuleCall_2_0_0; }
		
		//value=NumberValue
		public Assignment getValueAssignment_2_1() { return cValueAssignment_2_1; }
		
		//NumberValue
		public RuleCall getValueNumberValueParserRuleCall_2_1_0() { return cValueNumberValueParserRuleCall_2_1_0; }
		
		//value=BooleanValue
		public Assignment getValueAssignment_2_2() { return cValueAssignment_2_2; }
		
		//BooleanValue
		public RuleCall getValueBooleanValueParserRuleCall_2_2_0() { return cValueBooleanValueParserRuleCall_2_2_0; }
	}
	public class PropertyKeyElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.PropertyKey");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final RuleCall cSTRINGTerminalRuleCall_0 = (RuleCall)cAlternatives.eContents().get(0);
		private final RuleCall cIDTerminalRuleCall_1 = (RuleCall)cAlternatives.eContents().get(1);
		
		//PropertyKey IProperty hidden():
		//	STRING | ID;
		@Override public ParserRule getRule() { return rule; }
		
		//STRING | ID
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//STRING
		public RuleCall getSTRINGTerminalRuleCall_0() { return cSTRINGTerminalRuleCall_0; }
		
		//ID
		public RuleCall getIDTerminalRuleCall_1() { return cIDTerminalRuleCall_1; }
	}
	public class StringValueElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.StringValue");
		private final RuleCall cSTRINGTerminalRuleCall = (RuleCall)rule.eContents().get(1);
		
		//StringValue ecore::EJavaObject:
		//	STRING;
		@Override public ParserRule getRule() { return rule; }
		
		//STRING
		public RuleCall getSTRINGTerminalRuleCall() { return cSTRINGTerminalRuleCall; }
	}
	public class NumberValueElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.NumberValue");
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
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.BooleanValue");
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
	public class NumberElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.Number");
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
	public class JsonObjectElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.JsonObject");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Keyword cLeftCurlyBracketKeyword_0 = (Keyword)cGroup.eContents().get(0);
		private final Group cGroup_1 = (Group)cGroup.eContents().get(1);
		private final RuleCall cJsonMemberParserRuleCall_1_0 = (RuleCall)cGroup_1.eContents().get(0);
		private final Group cGroup_1_1 = (Group)cGroup_1.eContents().get(1);
		private final Keyword cCommaKeyword_1_1_0 = (Keyword)cGroup_1_1.eContents().get(0);
		private final RuleCall cJsonMemberParserRuleCall_1_1_1 = (RuleCall)cGroup_1_1.eContents().get(1);
		private final Keyword cCommaKeyword_2 = (Keyword)cGroup.eContents().get(2);
		private final Keyword cRightCurlyBracketKeyword_3 = (Keyword)cGroup.eContents().get(3);
		
		///* - - Arbitrary JSON (see https://www.json.org/json-en.html) - - */ JsonObject:
		//	'{' (JsonMember (',' JsonMember)*)? ','? '}';
		@Override public ParserRule getRule() { return rule; }
		
		//'{' (JsonMember (',' JsonMember)*)? ','? '}'
		public Group getGroup() { return cGroup; }
		
		//'{'
		public Keyword getLeftCurlyBracketKeyword_0() { return cLeftCurlyBracketKeyword_0; }
		
		//(JsonMember (',' JsonMember)*)?
		public Group getGroup_1() { return cGroup_1; }
		
		//JsonMember
		public RuleCall getJsonMemberParserRuleCall_1_0() { return cJsonMemberParserRuleCall_1_0; }
		
		//(',' JsonMember)*
		public Group getGroup_1_1() { return cGroup_1_1; }
		
		//','
		public Keyword getCommaKeyword_1_1_0() { return cCommaKeyword_1_1_0; }
		
		//JsonMember
		public RuleCall getJsonMemberParserRuleCall_1_1_1() { return cJsonMemberParserRuleCall_1_1_1; }
		
		//','?
		public Keyword getCommaKeyword_2() { return cCommaKeyword_2; }
		
		//'}'
		public Keyword getRightCurlyBracketKeyword_3() { return cRightCurlyBracketKeyword_3; }
	}
	public class JsonArrayElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.JsonArray");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Keyword cLeftSquareBracketKeyword_0 = (Keyword)cGroup.eContents().get(0);
		private final Group cGroup_1 = (Group)cGroup.eContents().get(1);
		private final RuleCall cJsonValueParserRuleCall_1_0 = (RuleCall)cGroup_1.eContents().get(0);
		private final Group cGroup_1_1 = (Group)cGroup_1.eContents().get(1);
		private final Keyword cCommaKeyword_1_1_0 = (Keyword)cGroup_1_1.eContents().get(0);
		private final RuleCall cJsonValueParserRuleCall_1_1_1 = (RuleCall)cGroup_1_1.eContents().get(1);
		private final Keyword cCommaKeyword_2 = (Keyword)cGroup.eContents().get(2);
		private final Keyword cRightSquareBracketKeyword_3 = (Keyword)cGroup.eContents().get(3);
		
		//JsonArray:
		//	'[' (JsonValue (',' JsonValue)*)? ','? ']';
		@Override public ParserRule getRule() { return rule; }
		
		//'[' (JsonValue (',' JsonValue)*)? ','? ']'
		public Group getGroup() { return cGroup; }
		
		//'['
		public Keyword getLeftSquareBracketKeyword_0() { return cLeftSquareBracketKeyword_0; }
		
		//(JsonValue (',' JsonValue)*)?
		public Group getGroup_1() { return cGroup_1; }
		
		//JsonValue
		public RuleCall getJsonValueParserRuleCall_1_0() { return cJsonValueParserRuleCall_1_0; }
		
		//(',' JsonValue)*
		public Group getGroup_1_1() { return cGroup_1_1; }
		
		//','
		public Keyword getCommaKeyword_1_1_0() { return cCommaKeyword_1_1_0; }
		
		//JsonValue
		public RuleCall getJsonValueParserRuleCall_1_1_1() { return cJsonValueParserRuleCall_1_1_1; }
		
		//','?
		public Keyword getCommaKeyword_2() { return cCommaKeyword_2; }
		
		//']'
		public Keyword getRightSquareBracketKeyword_3() { return cRightSquareBracketKeyword_3; }
	}
	public class JsonMemberElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.JsonMember");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Alternatives cAlternatives_0 = (Alternatives)cGroup.eContents().get(0);
		private final RuleCall cSTRINGTerminalRuleCall_0_0 = (RuleCall)cAlternatives_0.eContents().get(0);
		private final RuleCall cIDTerminalRuleCall_0_1 = (RuleCall)cAlternatives_0.eContents().get(1);
		private final Keyword cColonKeyword_1 = (Keyword)cGroup.eContents().get(1);
		private final RuleCall cJsonValueParserRuleCall_2 = (RuleCall)cGroup.eContents().get(2);
		
		//JsonMember:
		//	(STRING | ID) ':' JsonValue;
		@Override public ParserRule getRule() { return rule; }
		
		//(STRING | ID) ':' JsonValue
		public Group getGroup() { return cGroup; }
		
		//(STRING | ID)
		public Alternatives getAlternatives_0() { return cAlternatives_0; }
		
		//STRING
		public RuleCall getSTRINGTerminalRuleCall_0_0() { return cSTRINGTerminalRuleCall_0_0; }
		
		//ID
		public RuleCall getIDTerminalRuleCall_0_1() { return cIDTerminalRuleCall_0_1; }
		
		//':'
		public Keyword getColonKeyword_1() { return cColonKeyword_1; }
		
		//JsonValue
		public RuleCall getJsonValueParserRuleCall_2() { return cJsonValueParserRuleCall_2; }
	}
	public class JsonValueElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.JsonValue");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final RuleCall cSTRINGTerminalRuleCall_0 = (RuleCall)cAlternatives.eContents().get(0);
		private final RuleCall cNumberParserRuleCall_1 = (RuleCall)cAlternatives.eContents().get(1);
		private final RuleCall cJsonObjectParserRuleCall_2 = (RuleCall)cAlternatives.eContents().get(2);
		private final RuleCall cJsonArrayParserRuleCall_3 = (RuleCall)cAlternatives.eContents().get(3);
		private final RuleCall cBooleanValueParserRuleCall_4 = (RuleCall)cAlternatives.eContents().get(4);
		private final Keyword cNullKeyword_5 = (Keyword)cAlternatives.eContents().get(5);
		
		//JsonValue:
		//	STRING | Number | JsonObject | JsonArray | BooleanValue | 'null';
		@Override public ParserRule getRule() { return rule; }
		
		//STRING | Number | JsonObject | JsonArray | BooleanValue | 'null'
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//STRING
		public RuleCall getSTRINGTerminalRuleCall_0() { return cSTRINGTerminalRuleCall_0; }
		
		//Number
		public RuleCall getNumberParserRuleCall_1() { return cNumberParserRuleCall_1; }
		
		//JsonObject
		public RuleCall getJsonObjectParserRuleCall_2() { return cJsonObjectParserRuleCall_2; }
		
		//JsonArray
		public RuleCall getJsonArrayParserRuleCall_3() { return cJsonArrayParserRuleCall_3; }
		
		//BooleanValue
		public RuleCall getBooleanValueParserRuleCall_4() { return cBooleanValueParserRuleCall_4; }
		
		//'null'
		public Keyword getNullKeyword_5() { return cNullKeyword_5; }
	}
	public class KeyChildrenElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.KeyChildren");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final Keyword cChildrenKeyword_0 = (Keyword)cAlternatives.eContents().get(0);
		private final Keyword cChildrenKeyword_1 = (Keyword)cAlternatives.eContents().get(1);
		private final Keyword cChildrenKeyword_2 = (Keyword)cAlternatives.eContents().get(2);
		
		///* - - Keys (allowing single and double quotes) - - */ KeyChildren:
		//	'"children"' | "'children'" | 'children';
		@Override public ParserRule getRule() { return rule; }
		
		//'"children"' | "'children'" | 'children'
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//'"children"'
		public Keyword getChildrenKeyword_0() { return cChildrenKeyword_0; }
		
		//"'children'"
		public Keyword getChildrenKeyword_1() { return cChildrenKeyword_1; }
		
		//'children'
		public Keyword getChildrenKeyword_2() { return cChildrenKeyword_2; }
	}
	public class KeyPortsElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.KeyPorts");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final Keyword cPortsKeyword_0 = (Keyword)cAlternatives.eContents().get(0);
		private final Keyword cPortsKeyword_1 = (Keyword)cAlternatives.eContents().get(1);
		private final Keyword cPortsKeyword_2 = (Keyword)cAlternatives.eContents().get(2);
		
		//KeyPorts:
		//	'"ports"' | "'ports'" | 'ports';
		@Override public ParserRule getRule() { return rule; }
		
		//'"ports"' | "'ports'" | 'ports'
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//'"ports"'
		public Keyword getPortsKeyword_0() { return cPortsKeyword_0; }
		
		//"'ports'"
		public Keyword getPortsKeyword_1() { return cPortsKeyword_1; }
		
		//'ports'
		public Keyword getPortsKeyword_2() { return cPortsKeyword_2; }
	}
	public class KeyLabelsElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.KeyLabels");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final Keyword cLabelsKeyword_0 = (Keyword)cAlternatives.eContents().get(0);
		private final Keyword cLabelsKeyword_1 = (Keyword)cAlternatives.eContents().get(1);
		private final Keyword cLabelsKeyword_2 = (Keyword)cAlternatives.eContents().get(2);
		
		//KeyLabels:
		//	'"labels"' | "'labels'" | 'labels';
		@Override public ParserRule getRule() { return rule; }
		
		//'"labels"' | "'labels'" | 'labels'
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//'"labels"'
		public Keyword getLabelsKeyword_0() { return cLabelsKeyword_0; }
		
		//"'labels'"
		public Keyword getLabelsKeyword_1() { return cLabelsKeyword_1; }
		
		//'labels'
		public Keyword getLabelsKeyword_2() { return cLabelsKeyword_2; }
	}
	public class KeyEdgesElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.KeyEdges");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final Keyword cEdgesKeyword_0 = (Keyword)cAlternatives.eContents().get(0);
		private final Keyword cEdgesKeyword_1 = (Keyword)cAlternatives.eContents().get(1);
		private final Keyword cEdgesKeyword_2 = (Keyword)cAlternatives.eContents().get(2);
		
		//KeyEdges:
		//	'"edges"' | "'edges'" | 'edges';
		@Override public ParserRule getRule() { return rule; }
		
		//'"edges"' | "'edges'" | 'edges'
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//'"edges"'
		public Keyword getEdgesKeyword_0() { return cEdgesKeyword_0; }
		
		//"'edges'"
		public Keyword getEdgesKeyword_1() { return cEdgesKeyword_1; }
		
		//'edges'
		public Keyword getEdgesKeyword_2() { return cEdgesKeyword_2; }
	}
	public class KeyLayoutOptionsElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.KeyLayoutOptions");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final Keyword cLayoutOptionsKeyword_0 = (Keyword)cAlternatives.eContents().get(0);
		private final Keyword cLayoutOptionsKeyword_1 = (Keyword)cAlternatives.eContents().get(1);
		private final Keyword cLayoutOptionsKeyword_2 = (Keyword)cAlternatives.eContents().get(2);
		private final Keyword cPropertiesKeyword_3 = (Keyword)cAlternatives.eContents().get(3);
		private final Keyword cPropertiesKeyword_4 = (Keyword)cAlternatives.eContents().get(4);
		private final Keyword cPropertiesKeyword_5 = (Keyword)cAlternatives.eContents().get(5);
		
		//KeyLayoutOptions:
		//	'"layoutOptions"' | "'layoutOptions'" | 'layoutOptions'
		//	| '"properties"' | "'properties'" | 'properties';
		@Override public ParserRule getRule() { return rule; }
		
		//'"layoutOptions"' | "'layoutOptions'" | 'layoutOptions' | '"properties"' | "'properties'" | 'properties'
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//'"layoutOptions"'
		public Keyword getLayoutOptionsKeyword_0() { return cLayoutOptionsKeyword_0; }
		
		//"'layoutOptions'"
		public Keyword getLayoutOptionsKeyword_1() { return cLayoutOptionsKeyword_1; }
		
		//'layoutOptions'
		public Keyword getLayoutOptionsKeyword_2() { return cLayoutOptionsKeyword_2; }
		
		//'"properties"'
		public Keyword getPropertiesKeyword_3() { return cPropertiesKeyword_3; }
		
		//"'properties'"
		public Keyword getPropertiesKeyword_4() { return cPropertiesKeyword_4; }
		
		//'properties'
		public Keyword getPropertiesKeyword_5() { return cPropertiesKeyword_5; }
	}
	public class KeyIdElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.KeyId");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final Keyword cIdKeyword_0 = (Keyword)cAlternatives.eContents().get(0);
		private final Keyword cIdKeyword_1 = (Keyword)cAlternatives.eContents().get(1);
		private final Keyword cIdKeyword_2 = (Keyword)cAlternatives.eContents().get(2);
		
		//KeyId:
		//	'"id"' | "'id'" | 'id';
		@Override public ParserRule getRule() { return rule; }
		
		//'"id"' | "'id'" | 'id'
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//'"id"'
		public Keyword getIdKeyword_0() { return cIdKeyword_0; }
		
		//"'id'"
		public Keyword getIdKeyword_1() { return cIdKeyword_1; }
		
		//'id'
		public Keyword getIdKeyword_2() { return cIdKeyword_2; }
	}
	public class KeyXElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.KeyX");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final Keyword cXKeyword_0 = (Keyword)cAlternatives.eContents().get(0);
		private final Keyword cXKeyword_1 = (Keyword)cAlternatives.eContents().get(1);
		private final Keyword cXKeyword_2 = (Keyword)cAlternatives.eContents().get(2);
		
		//KeyX:
		//	'"x"' | "'x'" | 'x';
		@Override public ParserRule getRule() { return rule; }
		
		//'"x"' | "'x'" | 'x'
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//'"x"'
		public Keyword getXKeyword_0() { return cXKeyword_0; }
		
		//"'x'"
		public Keyword getXKeyword_1() { return cXKeyword_1; }
		
		//'x'
		public Keyword getXKeyword_2() { return cXKeyword_2; }
	}
	public class KeyYElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.KeyY");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final Keyword cYKeyword_0 = (Keyword)cAlternatives.eContents().get(0);
		private final Keyword cYKeyword_1 = (Keyword)cAlternatives.eContents().get(1);
		private final Keyword cYKeyword_2 = (Keyword)cAlternatives.eContents().get(2);
		
		//KeyY:
		//	'"y"' | "'y'" | 'y';
		@Override public ParserRule getRule() { return rule; }
		
		//'"y"' | "'y'" | 'y'
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//'"y"'
		public Keyword getYKeyword_0() { return cYKeyword_0; }
		
		//"'y'"
		public Keyword getYKeyword_1() { return cYKeyword_1; }
		
		//'y'
		public Keyword getYKeyword_2() { return cYKeyword_2; }
	}
	public class KeyWidthElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.KeyWidth");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final Keyword cWidthKeyword_0 = (Keyword)cAlternatives.eContents().get(0);
		private final Keyword cWidthKeyword_1 = (Keyword)cAlternatives.eContents().get(1);
		private final Keyword cWidthKeyword_2 = (Keyword)cAlternatives.eContents().get(2);
		
		//KeyWidth:
		//	'"width"' | "'width'" | 'width';
		@Override public ParserRule getRule() { return rule; }
		
		//'"width"' | "'width'" | 'width'
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//'"width"'
		public Keyword getWidthKeyword_0() { return cWidthKeyword_0; }
		
		//"'width'"
		public Keyword getWidthKeyword_1() { return cWidthKeyword_1; }
		
		//'width'
		public Keyword getWidthKeyword_2() { return cWidthKeyword_2; }
	}
	public class KeyHeightElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.KeyHeight");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final Keyword cHeightKeyword_0 = (Keyword)cAlternatives.eContents().get(0);
		private final Keyword cHeightKeyword_1 = (Keyword)cAlternatives.eContents().get(1);
		private final Keyword cHeightKeyword_2 = (Keyword)cAlternatives.eContents().get(2);
		
		//KeyHeight:
		//	'"height"' | "'height'" | 'height';
		@Override public ParserRule getRule() { return rule; }
		
		//'"height"' | "'height'" | 'height'
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//'"height"'
		public Keyword getHeightKeyword_0() { return cHeightKeyword_0; }
		
		//"'height'"
		public Keyword getHeightKeyword_1() { return cHeightKeyword_1; }
		
		//'height'
		public Keyword getHeightKeyword_2() { return cHeightKeyword_2; }
	}
	public class KeySourcesElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.KeySources");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final Keyword cSourcesKeyword_0 = (Keyword)cAlternatives.eContents().get(0);
		private final Keyword cSourcesKeyword_1 = (Keyword)cAlternatives.eContents().get(1);
		private final Keyword cSourcesKeyword_2 = (Keyword)cAlternatives.eContents().get(2);
		
		//KeySources:
		//	'"sources"' | "'sources'" | 'sources';
		@Override public ParserRule getRule() { return rule; }
		
		//'"sources"' | "'sources'" | 'sources'
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//'"sources"'
		public Keyword getSourcesKeyword_0() { return cSourcesKeyword_0; }
		
		//"'sources'"
		public Keyword getSourcesKeyword_1() { return cSourcesKeyword_1; }
		
		//'sources'
		public Keyword getSourcesKeyword_2() { return cSourcesKeyword_2; }
	}
	public class KeyTargetsElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.KeyTargets");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final Keyword cTargetsKeyword_0 = (Keyword)cAlternatives.eContents().get(0);
		private final Keyword cTargetsKeyword_1 = (Keyword)cAlternatives.eContents().get(1);
		private final Keyword cTargetsKeyword_2 = (Keyword)cAlternatives.eContents().get(2);
		
		//KeyTargets:
		//	'"targets"' | "'targets'" | 'targets';
		@Override public ParserRule getRule() { return rule; }
		
		//'"targets"' | "'targets'" | 'targets'
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//'"targets"'
		public Keyword getTargetsKeyword_0() { return cTargetsKeyword_0; }
		
		//"'targets'"
		public Keyword getTargetsKeyword_1() { return cTargetsKeyword_1; }
		
		//'targets'
		public Keyword getTargetsKeyword_2() { return cTargetsKeyword_2; }
	}
	public class KeyTextElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.KeyText");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final Keyword cTextKeyword_0 = (Keyword)cAlternatives.eContents().get(0);
		private final Keyword cTextKeyword_1 = (Keyword)cAlternatives.eContents().get(1);
		private final Keyword cTextKeyword_2 = (Keyword)cAlternatives.eContents().get(2);
		
		//KeyText:
		//	'"text"' | "'text'" | 'text';
		@Override public ParserRule getRule() { return rule; }
		
		//'"text"' | "'text'" | 'text'
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//'"text"'
		public Keyword getTextKeyword_0() { return cTextKeyword_0; }
		
		//"'text'"
		public Keyword getTextKeyword_1() { return cTextKeyword_1; }
		
		//'text'
		public Keyword getTextKeyword_2() { return cTextKeyword_2; }
	}
	
	
	private final ElkNodeElements pElkNode;
	private final NodeElementElements pNodeElement;
	private final ElkPortElements pElkPort;
	private final PortElementElements pPortElement;
	private final ElkLabelElements pElkLabel;
	private final LabelElementElements pLabelElement;
	private final ElkEdgeElements pElkEdge;
	private final EdgeElementElements pEdgeElement;
	private final ElkEdgeSourcesElements pElkEdgeSources;
	private final ElkEdgeTargetsElements pElkEdgeTargets;
	private final ElkIdElements pElkId;
	private final ElkNodeChildrenElements pElkNodeChildren;
	private final ElkNodePortsElements pElkNodePorts;
	private final ElkNodeEdgesElements pElkNodeEdges;
	private final ElkGraphElementLabelsElements pElkGraphElementLabels;
	private final ElkGraphElementPropertiesElements pElkGraphElementProperties;
	private final ShapeElementElements pShapeElement;
	private final PropertyElements pProperty;
	private final PropertyKeyElements pPropertyKey;
	private final StringValueElements pStringValue;
	private final NumberValueElements pNumberValue;
	private final BooleanValueElements pBooleanValue;
	private final NumberElements pNumber;
	private final TerminalRule tSIGNED_INT;
	private final TerminalRule tFLOAT;
	private final JsonObjectElements pJsonObject;
	private final JsonArrayElements pJsonArray;
	private final JsonMemberElements pJsonMember;
	private final JsonValueElements pJsonValue;
	private final KeyChildrenElements pKeyChildren;
	private final KeyPortsElements pKeyPorts;
	private final KeyLabelsElements pKeyLabels;
	private final KeyEdgesElements pKeyEdges;
	private final KeyLayoutOptionsElements pKeyLayoutOptions;
	private final KeyIdElements pKeyId;
	private final KeyXElements pKeyX;
	private final KeyYElements pKeyY;
	private final KeyWidthElements pKeyWidth;
	private final KeyHeightElements pKeyHeight;
	private final KeySourcesElements pKeySources;
	private final KeyTargetsElements pKeyTargets;
	private final KeyTextElements pKeyText;
	
	private final Grammar grammar;
	
	private final TerminalsGrammarAccess gaTerminals;

	@Inject
	public ElkGraphJsonGrammarAccess(GrammarProvider grammarProvider,
			TerminalsGrammarAccess gaTerminals) {
		this.grammar = internalFindGrammar(grammarProvider);
		this.gaTerminals = gaTerminals;
		this.pElkNode = new ElkNodeElements();
		this.pNodeElement = new NodeElementElements();
		this.pElkPort = new ElkPortElements();
		this.pPortElement = new PortElementElements();
		this.pElkLabel = new ElkLabelElements();
		this.pLabelElement = new LabelElementElements();
		this.pElkEdge = new ElkEdgeElements();
		this.pEdgeElement = new EdgeElementElements();
		this.pElkEdgeSources = new ElkEdgeSourcesElements();
		this.pElkEdgeTargets = new ElkEdgeTargetsElements();
		this.pElkId = new ElkIdElements();
		this.pElkNodeChildren = new ElkNodeChildrenElements();
		this.pElkNodePorts = new ElkNodePortsElements();
		this.pElkNodeEdges = new ElkNodeEdgesElements();
		this.pElkGraphElementLabels = new ElkGraphElementLabelsElements();
		this.pElkGraphElementProperties = new ElkGraphElementPropertiesElements();
		this.pShapeElement = new ShapeElementElements();
		this.pProperty = new PropertyElements();
		this.pPropertyKey = new PropertyKeyElements();
		this.pStringValue = new StringValueElements();
		this.pNumberValue = new NumberValueElements();
		this.pBooleanValue = new BooleanValueElements();
		this.pNumber = new NumberElements();
		this.tSIGNED_INT = (TerminalRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.SIGNED_INT");
		this.tFLOAT = (TerminalRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.graph.json.text.ElkGraphJson.FLOAT");
		this.pJsonObject = new JsonObjectElements();
		this.pJsonArray = new JsonArrayElements();
		this.pJsonMember = new JsonMemberElements();
		this.pJsonValue = new JsonValueElements();
		this.pKeyChildren = new KeyChildrenElements();
		this.pKeyPorts = new KeyPortsElements();
		this.pKeyLabels = new KeyLabelsElements();
		this.pKeyEdges = new KeyEdgesElements();
		this.pKeyLayoutOptions = new KeyLayoutOptionsElements();
		this.pKeyId = new KeyIdElements();
		this.pKeyX = new KeyXElements();
		this.pKeyY = new KeyYElements();
		this.pKeyWidth = new KeyWidthElements();
		this.pKeyHeight = new KeyHeightElements();
		this.pKeySources = new KeySourcesElements();
		this.pKeyTargets = new KeyTargetsElements();
		this.pKeyText = new KeyTextElements();
	}
	
	protected Grammar internalFindGrammar(GrammarProvider grammarProvider) {
		Grammar grammar = grammarProvider.getGrammar(this);
		while (grammar != null) {
			if ("org.eclipse.elk.graph.json.text.ElkGraphJson".equals(grammar.getName())) {
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

	
	///* - - Nodes - - */ ElkNode:
	//	{ElkNode} '{' (NodeElement (',' NodeElement)*)? ','? '}';
	public ElkNodeElements getElkNodeAccess() {
		return pElkNode;
	}
	
	public ParserRule getElkNodeRule() {
		return getElkNodeAccess().getRule();
	}
	
	//fragment NodeElement returns ElkNode:
	//	ElkId
	//	| ShapeElement
	//	| KeyChildren ':' ElkNodeChildren
	//	| KeyPorts ':' ElkNodePorts
	//	| KeyLabels ':' ElkGraphElementLabels
	//	| KeyEdges ':' ElkNodeEdges
	//	| KeyLayoutOptions ':' ElkGraphElementProperties
	//	| JsonMember;
	public NodeElementElements getNodeElementAccess() {
		return pNodeElement;
	}
	
	public ParserRule getNodeElementRule() {
		return getNodeElementAccess().getRule();
	}
	
	///* - - Ports - - */ ElkPort:
	//	'{' PortElement (',' PortElement)* ','? '}';
	public ElkPortElements getElkPortAccess() {
		return pElkPort;
	}
	
	public ParserRule getElkPortRule() {
		return getElkPortAccess().getRule();
	}
	
	//fragment PortElement returns ElkPort:
	//	ElkId
	//	| ShapeElement
	//	| KeyLabels ':' ElkGraphElementLabels
	//	| KeyLayoutOptions ':' ElkGraphElementProperties
	//	| JsonMember;
	public PortElementElements getPortElementAccess() {
		return pPortElement;
	}
	
	public ParserRule getPortElementRule() {
		return getPortElementAccess().getRule();
	}
	
	///* - - Labels - - */ ElkLabel:
	//	'{' LabelElement (',' LabelElement)* ','? '}';
	public ElkLabelElements getElkLabelAccess() {
		return pElkLabel;
	}
	
	public ParserRule getElkLabelRule() {
		return getElkLabelAccess().getRule();
	}
	
	//fragment LabelElement returns ElkLabel:
	//	ElkId
	//	| ShapeElement
	//	| KeyText ':' text=STRING
	//	| KeyLabels ':' ElkGraphElementLabels
	//	| KeyLayoutOptions ':' ElkGraphElementProperties
	//	| JsonMember;
	public LabelElementElements getLabelElementAccess() {
		return pLabelElement;
	}
	
	public ParserRule getLabelElementRule() {
		return getLabelElementAccess().getRule();
	}
	
	///* - - Edges - - */ ElkEdge:
	//	'{' EdgeElement (',' EdgeElement)* ','? '}';
	public ElkEdgeElements getElkEdgeAccess() {
		return pElkEdge;
	}
	
	public ParserRule getElkEdgeRule() {
		return getElkEdgeAccess().getRule();
	}
	
	//fragment EdgeElement returns ElkEdge:
	//	ElkId
	//	// TODO EdgeLayout
	//	| KeySources ':' ElkEdgeSources
	//	| KeyTargets ':' ElkEdgeTargets
	//	| KeyLabels ':' ElkGraphElementLabels
	//	| KeyLayoutOptions ':' ElkGraphElementProperties
	//	| JsonMember;
	public EdgeElementElements getEdgeElementAccess() {
		return pEdgeElement;
	}
	
	public ParserRule getEdgeElementRule() {
		return getEdgeElementAccess().getRule();
	}
	
	///* SuppressWarnings[BidirectionalReference] */ fragment ElkEdgeSources returns ElkEdge:
	//	'[' (sources+=[ElkConnectableShape|STRING] (',' sources+=[ElkConnectableShape|STRING])*)? ','? ']';
	public ElkEdgeSourcesElements getElkEdgeSourcesAccess() {
		return pElkEdgeSources;
	}
	
	public ParserRule getElkEdgeSourcesRule() {
		return getElkEdgeSourcesAccess().getRule();
	}
	
	///* SuppressWarnings[BidirectionalReference] */ fragment ElkEdgeTargets returns ElkEdge:
	//	'[' (targets+=[ElkConnectableShape|STRING] (',' targets+=[ElkConnectableShape|STRING])*)? ','? ']';
	public ElkEdgeTargetsElements getElkEdgeTargetsAccess() {
		return pElkEdgeTargets;
	}
	
	public ParserRule getElkEdgeTargetsRule() {
		return getElkEdgeTargetsAccess().getRule();
	}
	
	///* - - ElkGraphElement fragments - - */ fragment ElkId returns ElkGraphElement:
	//	KeyId ':' identifier=STRING;
	public ElkIdElements getElkIdAccess() {
		return pElkId;
	}
	
	public ParserRule getElkIdRule() {
		return getElkIdAccess().getRule();
	}
	
	//fragment ElkNodeChildren returns ElkNode:
	//	'[' (children+=ElkNode (',' children+=ElkNode)*)? ','? ']';
	public ElkNodeChildrenElements getElkNodeChildrenAccess() {
		return pElkNodeChildren;
	}
	
	public ParserRule getElkNodeChildrenRule() {
		return getElkNodeChildrenAccess().getRule();
	}
	
	//fragment ElkNodePorts returns ElkNode:
	//	'[' (ports+=ElkPort (',' ports+=ElkPort)*)? ','? ']';
	public ElkNodePortsElements getElkNodePortsAccess() {
		return pElkNodePorts;
	}
	
	public ParserRule getElkNodePortsRule() {
		return getElkNodePortsAccess().getRule();
	}
	
	//fragment ElkNodeEdges returns ElkNode:
	//	'[' (containedEdges+=ElkEdge (',' containedEdges+=ElkEdge)*)? ','? ']';
	public ElkNodeEdgesElements getElkNodeEdgesAccess() {
		return pElkNodeEdges;
	}
	
	public ParserRule getElkNodeEdgesRule() {
		return getElkNodeEdgesAccess().getRule();
	}
	
	//fragment ElkGraphElementLabels returns ElkGraphElement:
	//	'[' (labels+=ElkLabel (',' labels+=ElkLabel)*)? ','? ']';
	public ElkGraphElementLabelsElements getElkGraphElementLabelsAccess() {
		return pElkGraphElementLabels;
	}
	
	public ParserRule getElkGraphElementLabelsRule() {
		return getElkGraphElementLabelsAccess().getRule();
	}
	
	//fragment ElkGraphElementProperties returns ElkGraphElement:
	//	'{' (properties+=Property (',' properties+=Property)*)? ','? '}';
	public ElkGraphElementPropertiesElements getElkGraphElementPropertiesAccess() {
		return pElkGraphElementProperties;
	}
	
	public ParserRule getElkGraphElementPropertiesRule() {
		return getElkGraphElementPropertiesAccess().getRule();
	}
	
	///* - - Shape - - */ fragment ShapeElement returns ElkShape:
	//	KeyX ':' x=Number
	//	| KeyY ':' y=Number
	//	| KeyWidth ':' width=Number
	//	| KeyHeight ':' height=Number;
	public ShapeElementElements getShapeElementAccess() {
		return pShapeElement;
	}
	
	public ParserRule getShapeElementRule() {
		return getShapeElementAccess().getRule();
	}
	
	//Property ElkPropertyToValueMapEntry:
	//	key=PropertyKey ':' (value=StringValue | value=NumberValue | value=BooleanValue);
	public PropertyElements getPropertyAccess() {
		return pProperty;
	}
	
	public ParserRule getPropertyRule() {
		return getPropertyAccess().getRule();
	}
	
	//PropertyKey IProperty hidden():
	//	STRING | ID;
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
	
	//Number ecore::EDouble:
	//	SIGNED_INT | FLOAT;
	public NumberElements getNumberAccess() {
		return pNumber;
	}
	
	public ParserRule getNumberRule() {
		return getNumberAccess().getRule();
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
	
	///* - - Arbitrary JSON (see https://www.json.org/json-en.html) - - */ JsonObject:
	//	'{' (JsonMember (',' JsonMember)*)? ','? '}';
	public JsonObjectElements getJsonObjectAccess() {
		return pJsonObject;
	}
	
	public ParserRule getJsonObjectRule() {
		return getJsonObjectAccess().getRule();
	}
	
	//JsonArray:
	//	'[' (JsonValue (',' JsonValue)*)? ','? ']';
	public JsonArrayElements getJsonArrayAccess() {
		return pJsonArray;
	}
	
	public ParserRule getJsonArrayRule() {
		return getJsonArrayAccess().getRule();
	}
	
	//JsonMember:
	//	(STRING | ID) ':' JsonValue;
	public JsonMemberElements getJsonMemberAccess() {
		return pJsonMember;
	}
	
	public ParserRule getJsonMemberRule() {
		return getJsonMemberAccess().getRule();
	}
	
	//JsonValue:
	//	STRING | Number | JsonObject | JsonArray | BooleanValue | 'null';
	public JsonValueElements getJsonValueAccess() {
		return pJsonValue;
	}
	
	public ParserRule getJsonValueRule() {
		return getJsonValueAccess().getRule();
	}
	
	///* - - Keys (allowing single and double quotes) - - */ KeyChildren:
	//	'"children"' | "'children'" | 'children';
	public KeyChildrenElements getKeyChildrenAccess() {
		return pKeyChildren;
	}
	
	public ParserRule getKeyChildrenRule() {
		return getKeyChildrenAccess().getRule();
	}
	
	//KeyPorts:
	//	'"ports"' | "'ports'" | 'ports';
	public KeyPortsElements getKeyPortsAccess() {
		return pKeyPorts;
	}
	
	public ParserRule getKeyPortsRule() {
		return getKeyPortsAccess().getRule();
	}
	
	//KeyLabels:
	//	'"labels"' | "'labels'" | 'labels';
	public KeyLabelsElements getKeyLabelsAccess() {
		return pKeyLabels;
	}
	
	public ParserRule getKeyLabelsRule() {
		return getKeyLabelsAccess().getRule();
	}
	
	//KeyEdges:
	//	'"edges"' | "'edges'" | 'edges';
	public KeyEdgesElements getKeyEdgesAccess() {
		return pKeyEdges;
	}
	
	public ParserRule getKeyEdgesRule() {
		return getKeyEdgesAccess().getRule();
	}
	
	//KeyLayoutOptions:
	//	'"layoutOptions"' | "'layoutOptions'" | 'layoutOptions'
	//	| '"properties"' | "'properties'" | 'properties';
	public KeyLayoutOptionsElements getKeyLayoutOptionsAccess() {
		return pKeyLayoutOptions;
	}
	
	public ParserRule getKeyLayoutOptionsRule() {
		return getKeyLayoutOptionsAccess().getRule();
	}
	
	//KeyId:
	//	'"id"' | "'id'" | 'id';
	public KeyIdElements getKeyIdAccess() {
		return pKeyId;
	}
	
	public ParserRule getKeyIdRule() {
		return getKeyIdAccess().getRule();
	}
	
	//KeyX:
	//	'"x"' | "'x'" | 'x';
	public KeyXElements getKeyXAccess() {
		return pKeyX;
	}
	
	public ParserRule getKeyXRule() {
		return getKeyXAccess().getRule();
	}
	
	//KeyY:
	//	'"y"' | "'y'" | 'y';
	public KeyYElements getKeyYAccess() {
		return pKeyY;
	}
	
	public ParserRule getKeyYRule() {
		return getKeyYAccess().getRule();
	}
	
	//KeyWidth:
	//	'"width"' | "'width'" | 'width';
	public KeyWidthElements getKeyWidthAccess() {
		return pKeyWidth;
	}
	
	public ParserRule getKeyWidthRule() {
		return getKeyWidthAccess().getRule();
	}
	
	//KeyHeight:
	//	'"height"' | "'height'" | 'height';
	public KeyHeightElements getKeyHeightAccess() {
		return pKeyHeight;
	}
	
	public ParserRule getKeyHeightRule() {
		return getKeyHeightAccess().getRule();
	}
	
	//KeySources:
	//	'"sources"' | "'sources'" | 'sources';
	public KeySourcesElements getKeySourcesAccess() {
		return pKeySources;
	}
	
	public ParserRule getKeySourcesRule() {
		return getKeySourcesAccess().getRule();
	}
	
	//KeyTargets:
	//	'"targets"' | "'targets'" | 'targets';
	public KeyTargetsElements getKeyTargetsAccess() {
		return pKeyTargets;
	}
	
	public ParserRule getKeyTargetsRule() {
		return getKeyTargetsAccess().getRule();
	}
	
	//KeyText:
	//	'"text"' | "'text'" | 'text';
	public KeyTextElements getKeyTextAccess() {
		return pKeyText;
	}
	
	public ParserRule getKeyTextRule() {
		return getKeyTextAccess().getRule();
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
	//	'"' ('\\' . | !('\\' | '"'))* '"' | "'" ('\\' . | !('\\' | "'"))* "'";
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
