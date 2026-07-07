/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 *
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
grammar InternalElkGraphJson;

options {
	superClass=AbstractInternalAntlrParser;
}

@lexer::header {
package org.eclipse.elk.graph.json.text.parser.antlr.internal;

// Hack: Use our own Lexer superclass by means of import. 
// Currently there is no other way to specify the superclass for the lexer.
import org.eclipse.xtext.parser.antlr.Lexer;
}

@parser::header {
package org.eclipse.elk.graph.json.text.parser.antlr.internal;

import org.eclipse.xtext.*;
import org.eclipse.xtext.parser.*;
import org.eclipse.xtext.parser.impl.*;
import org.eclipse.emf.ecore.util.EcoreUtil;
import org.eclipse.emf.ecore.EObject;
import org.eclipse.xtext.parser.antlr.AbstractInternalAntlrParser;
import org.eclipse.xtext.parser.antlr.XtextTokenStream;
import org.eclipse.xtext.parser.antlr.XtextTokenStream.HiddenTokens;
import org.eclipse.xtext.parser.antlr.AntlrDatatypeRuleToken;
import org.eclipse.elk.graph.json.text.services.ElkGraphJsonGrammarAccess;

}

@parser::members {

 	private ElkGraphJsonGrammarAccess grammarAccess;

    public InternalElkGraphJsonParser(TokenStream input, ElkGraphJsonGrammarAccess grammarAccess) {
        this(input);
        this.grammarAccess = grammarAccess;
        registerRules(grammarAccess.getGrammar());
    }

    @Override
    protected String getFirstRuleName() {
    	return "ElkNode";
   	}

   	@Override
   	protected ElkGraphJsonGrammarAccess getGrammarAccess() {
   		return grammarAccess;
   	}

}

@rulecatch {
    catch (RecognitionException re) {
        recover(input,re);
        appendSkippedTokens();
    }
}

// Entry rule entryRuleElkNode
entryRuleElkNode returns [EObject current=null]:
	{ newCompositeNode(grammarAccess.getElkNodeRule()); }
	iv_ruleElkNode=ruleElkNode
	{ $current=$iv_ruleElkNode.current; }
	EOF;

// Rule ElkNode
ruleElkNode returns [EObject current=null]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		(
			{
				$current = forceCreateModelElement(
					grammarAccess.getElkNodeAccess().getElkNodeAction_0(),
					$current);
			}
		)
		otherlv_1='{'
		{
			newLeafNode(otherlv_1, grammarAccess.getElkNodeAccess().getLeftCurlyBracketKeyword_1());
		}
		(
			{
				if ($current==null) {
					$current = createModelElement(grammarAccess.getElkNodeRule());
				}
				newCompositeNode(grammarAccess.getElkNodeAccess().getNodeElementParserRuleCall_2_0());
			}
			this_NodeElement_2=ruleNodeElement[$current]
			{
				$current = $this_NodeElement_2.current;
				afterParserOrEnumRuleCall();
			}
			(
				otherlv_3=','
				{
					newLeafNode(otherlv_3, grammarAccess.getElkNodeAccess().getCommaKeyword_2_1_0());
				}
				{
					if ($current==null) {
						$current = createModelElement(grammarAccess.getElkNodeRule());
					}
					newCompositeNode(grammarAccess.getElkNodeAccess().getNodeElementParserRuleCall_2_1_1());
				}
				this_NodeElement_4=ruleNodeElement[$current]
				{
					$current = $this_NodeElement_4.current;
					afterParserOrEnumRuleCall();
				}
			)*
		)?
		(
			otherlv_5=','
			{
				newLeafNode(otherlv_5, grammarAccess.getElkNodeAccess().getCommaKeyword_3());
			}
		)?
		otherlv_6='}'
		{
			newLeafNode(otherlv_6, grammarAccess.getElkNodeAccess().getRightCurlyBracketKeyword_4());
		}
	)
;


// Rule NodeElement
ruleNodeElement[EObject in_current]  returns [EObject current=in_current]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		{
			if ($current==null) {
				$current = createModelElement(grammarAccess.getNodeElementRule());
			}
			newCompositeNode(grammarAccess.getNodeElementAccess().getElkIdParserRuleCall_0());
		}
		this_ElkId_0=ruleElkId[$current]
		{
			$current = $this_ElkId_0.current;
			afterParserOrEnumRuleCall();
		}
		    |
		{
			if ($current==null) {
				$current = createModelElement(grammarAccess.getNodeElementRule());
			}
			newCompositeNode(grammarAccess.getNodeElementAccess().getShapeElementParserRuleCall_1());
		}
		this_ShapeElement_1=ruleShapeElement[$current]
		{
			$current = $this_ShapeElement_1.current;
			afterParserOrEnumRuleCall();
		}
		    |
		(
			{
				newCompositeNode(grammarAccess.getNodeElementAccess().getKeyChildrenParserRuleCall_2_0());
			}
			ruleKeyChildren
			{
				afterParserOrEnumRuleCall();
			}
			otherlv_3=':'
			{
				newLeafNode(otherlv_3, grammarAccess.getNodeElementAccess().getColonKeyword_2_1());
			}
			{
				if ($current==null) {
					$current = createModelElement(grammarAccess.getNodeElementRule());
				}
				newCompositeNode(grammarAccess.getNodeElementAccess().getElkNodeChildrenParserRuleCall_2_2());
			}
			this_ElkNodeChildren_4=ruleElkNodeChildren[$current]
			{
				$current = $this_ElkNodeChildren_4.current;
				afterParserOrEnumRuleCall();
			}
		)
		    |
		(
			{
				newCompositeNode(grammarAccess.getNodeElementAccess().getKeyPortsParserRuleCall_3_0());
			}
			ruleKeyPorts
			{
				afterParserOrEnumRuleCall();
			}
			otherlv_6=':'
			{
				newLeafNode(otherlv_6, grammarAccess.getNodeElementAccess().getColonKeyword_3_1());
			}
			{
				if ($current==null) {
					$current = createModelElement(grammarAccess.getNodeElementRule());
				}
				newCompositeNode(grammarAccess.getNodeElementAccess().getElkNodePortsParserRuleCall_3_2());
			}
			this_ElkNodePorts_7=ruleElkNodePorts[$current]
			{
				$current = $this_ElkNodePorts_7.current;
				afterParserOrEnumRuleCall();
			}
		)
		    |
		(
			{
				newCompositeNode(grammarAccess.getNodeElementAccess().getKeyLabelsParserRuleCall_4_0());
			}
			ruleKeyLabels
			{
				afterParserOrEnumRuleCall();
			}
			otherlv_9=':'
			{
				newLeafNode(otherlv_9, grammarAccess.getNodeElementAccess().getColonKeyword_4_1());
			}
			{
				if ($current==null) {
					$current = createModelElement(grammarAccess.getNodeElementRule());
				}
				newCompositeNode(grammarAccess.getNodeElementAccess().getElkGraphElementLabelsParserRuleCall_4_2());
			}
			this_ElkGraphElementLabels_10=ruleElkGraphElementLabels[$current]
			{
				$current = $this_ElkGraphElementLabels_10.current;
				afterParserOrEnumRuleCall();
			}
		)
		    |
		(
			{
				newCompositeNode(grammarAccess.getNodeElementAccess().getKeyEdgesParserRuleCall_5_0());
			}
			ruleKeyEdges
			{
				afterParserOrEnumRuleCall();
			}
			otherlv_12=':'
			{
				newLeafNode(otherlv_12, grammarAccess.getNodeElementAccess().getColonKeyword_5_1());
			}
			{
				if ($current==null) {
					$current = createModelElement(grammarAccess.getNodeElementRule());
				}
				newCompositeNode(grammarAccess.getNodeElementAccess().getElkNodeEdgesParserRuleCall_5_2());
			}
			this_ElkNodeEdges_13=ruleElkNodeEdges[$current]
			{
				$current = $this_ElkNodeEdges_13.current;
				afterParserOrEnumRuleCall();
			}
		)
		    |
		(
			{
				newCompositeNode(grammarAccess.getNodeElementAccess().getKeyLayoutOptionsParserRuleCall_6_0());
			}
			ruleKeyLayoutOptions
			{
				afterParserOrEnumRuleCall();
			}
			otherlv_15=':'
			{
				newLeafNode(otherlv_15, grammarAccess.getNodeElementAccess().getColonKeyword_6_1());
			}
			{
				if ($current==null) {
					$current = createModelElement(grammarAccess.getNodeElementRule());
				}
				newCompositeNode(grammarAccess.getNodeElementAccess().getElkGraphElementPropertiesParserRuleCall_6_2());
			}
			this_ElkGraphElementProperties_16=ruleElkGraphElementProperties[$current]
			{
				$current = $this_ElkGraphElementProperties_16.current;
				afterParserOrEnumRuleCall();
			}
		)
		    |
		{
			newCompositeNode(grammarAccess.getNodeElementAccess().getJsonMemberParserRuleCall_7());
		}
		ruleJsonMember
		{
			afterParserOrEnumRuleCall();
		}
	)
;

// Entry rule entryRuleElkPort
entryRuleElkPort returns [EObject current=null]:
	{ newCompositeNode(grammarAccess.getElkPortRule()); }
	iv_ruleElkPort=ruleElkPort
	{ $current=$iv_ruleElkPort.current; }
	EOF;

// Rule ElkPort
ruleElkPort returns [EObject current=null]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		otherlv_0='{'
		{
			newLeafNode(otherlv_0, grammarAccess.getElkPortAccess().getLeftCurlyBracketKeyword_0());
		}
		{
			if ($current==null) {
				$current = createModelElement(grammarAccess.getElkPortRule());
			}
			newCompositeNode(grammarAccess.getElkPortAccess().getPortElementParserRuleCall_1());
		}
		this_PortElement_1=rulePortElement[$current]
		{
			$current = $this_PortElement_1.current;
			afterParserOrEnumRuleCall();
		}
		(
			otherlv_2=','
			{
				newLeafNode(otherlv_2, grammarAccess.getElkPortAccess().getCommaKeyword_2_0());
			}
			{
				if ($current==null) {
					$current = createModelElement(grammarAccess.getElkPortRule());
				}
				newCompositeNode(grammarAccess.getElkPortAccess().getPortElementParserRuleCall_2_1());
			}
			this_PortElement_3=rulePortElement[$current]
			{
				$current = $this_PortElement_3.current;
				afterParserOrEnumRuleCall();
			}
		)*
		(
			otherlv_4=','
			{
				newLeafNode(otherlv_4, grammarAccess.getElkPortAccess().getCommaKeyword_3());
			}
		)?
		otherlv_5='}'
		{
			newLeafNode(otherlv_5, grammarAccess.getElkPortAccess().getRightCurlyBracketKeyword_4());
		}
	)
;


// Rule PortElement
rulePortElement[EObject in_current]  returns [EObject current=in_current]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		{
			if ($current==null) {
				$current = createModelElement(grammarAccess.getPortElementRule());
			}
			newCompositeNode(grammarAccess.getPortElementAccess().getElkIdParserRuleCall_0());
		}
		this_ElkId_0=ruleElkId[$current]
		{
			$current = $this_ElkId_0.current;
			afterParserOrEnumRuleCall();
		}
		    |
		{
			if ($current==null) {
				$current = createModelElement(grammarAccess.getPortElementRule());
			}
			newCompositeNode(grammarAccess.getPortElementAccess().getShapeElementParserRuleCall_1());
		}
		this_ShapeElement_1=ruleShapeElement[$current]
		{
			$current = $this_ShapeElement_1.current;
			afterParserOrEnumRuleCall();
		}
		    |
		(
			{
				newCompositeNode(grammarAccess.getPortElementAccess().getKeyLabelsParserRuleCall_2_0());
			}
			ruleKeyLabels
			{
				afterParserOrEnumRuleCall();
			}
			otherlv_3=':'
			{
				newLeafNode(otherlv_3, grammarAccess.getPortElementAccess().getColonKeyword_2_1());
			}
			{
				if ($current==null) {
					$current = createModelElement(grammarAccess.getPortElementRule());
				}
				newCompositeNode(grammarAccess.getPortElementAccess().getElkGraphElementLabelsParserRuleCall_2_2());
			}
			this_ElkGraphElementLabels_4=ruleElkGraphElementLabels[$current]
			{
				$current = $this_ElkGraphElementLabels_4.current;
				afterParserOrEnumRuleCall();
			}
		)
		    |
		(
			{
				newCompositeNode(grammarAccess.getPortElementAccess().getKeyLayoutOptionsParserRuleCall_3_0());
			}
			ruleKeyLayoutOptions
			{
				afterParserOrEnumRuleCall();
			}
			otherlv_6=':'
			{
				newLeafNode(otherlv_6, grammarAccess.getPortElementAccess().getColonKeyword_3_1());
			}
			{
				if ($current==null) {
					$current = createModelElement(grammarAccess.getPortElementRule());
				}
				newCompositeNode(grammarAccess.getPortElementAccess().getElkGraphElementPropertiesParserRuleCall_3_2());
			}
			this_ElkGraphElementProperties_7=ruleElkGraphElementProperties[$current]
			{
				$current = $this_ElkGraphElementProperties_7.current;
				afterParserOrEnumRuleCall();
			}
		)
		    |
		{
			newCompositeNode(grammarAccess.getPortElementAccess().getJsonMemberParserRuleCall_4());
		}
		ruleJsonMember
		{
			afterParserOrEnumRuleCall();
		}
	)
;

// Entry rule entryRuleElkLabel
entryRuleElkLabel returns [EObject current=null]:
	{ newCompositeNode(grammarAccess.getElkLabelRule()); }
	iv_ruleElkLabel=ruleElkLabel
	{ $current=$iv_ruleElkLabel.current; }
	EOF;

// Rule ElkLabel
ruleElkLabel returns [EObject current=null]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		otherlv_0='{'
		{
			newLeafNode(otherlv_0, grammarAccess.getElkLabelAccess().getLeftCurlyBracketKeyword_0());
		}
		{
			if ($current==null) {
				$current = createModelElement(grammarAccess.getElkLabelRule());
			}
			newCompositeNode(grammarAccess.getElkLabelAccess().getLabelElementParserRuleCall_1());
		}
		this_LabelElement_1=ruleLabelElement[$current]
		{
			$current = $this_LabelElement_1.current;
			afterParserOrEnumRuleCall();
		}
		(
			otherlv_2=','
			{
				newLeafNode(otherlv_2, grammarAccess.getElkLabelAccess().getCommaKeyword_2_0());
			}
			{
				if ($current==null) {
					$current = createModelElement(grammarAccess.getElkLabelRule());
				}
				newCompositeNode(grammarAccess.getElkLabelAccess().getLabelElementParserRuleCall_2_1());
			}
			this_LabelElement_3=ruleLabelElement[$current]
			{
				$current = $this_LabelElement_3.current;
				afterParserOrEnumRuleCall();
			}
		)*
		(
			otherlv_4=','
			{
				newLeafNode(otherlv_4, grammarAccess.getElkLabelAccess().getCommaKeyword_3());
			}
		)?
		otherlv_5='}'
		{
			newLeafNode(otherlv_5, grammarAccess.getElkLabelAccess().getRightCurlyBracketKeyword_4());
		}
	)
;


// Rule LabelElement
ruleLabelElement[EObject in_current]  returns [EObject current=in_current]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		{
			if ($current==null) {
				$current = createModelElement(grammarAccess.getLabelElementRule());
			}
			newCompositeNode(grammarAccess.getLabelElementAccess().getElkIdParserRuleCall_0());
		}
		this_ElkId_0=ruleElkId[$current]
		{
			$current = $this_ElkId_0.current;
			afterParserOrEnumRuleCall();
		}
		    |
		{
			if ($current==null) {
				$current = createModelElement(grammarAccess.getLabelElementRule());
			}
			newCompositeNode(grammarAccess.getLabelElementAccess().getShapeElementParserRuleCall_1());
		}
		this_ShapeElement_1=ruleShapeElement[$current]
		{
			$current = $this_ShapeElement_1.current;
			afterParserOrEnumRuleCall();
		}
		    |
		(
			{
				newCompositeNode(grammarAccess.getLabelElementAccess().getKeyTextParserRuleCall_2_0());
			}
			ruleKeyText
			{
				afterParserOrEnumRuleCall();
			}
			otherlv_3=':'
			{
				newLeafNode(otherlv_3, grammarAccess.getLabelElementAccess().getColonKeyword_2_1());
			}
			(
				(
					lv_text_4_0=RULE_STRING
					{
						newLeafNode(lv_text_4_0, grammarAccess.getLabelElementAccess().getTextSTRINGTerminalRuleCall_2_2_0());
					}
					{
						if ($current==null) {
							$current = createModelElement(grammarAccess.getLabelElementRule());
						}
						setWithLastConsumed(
							$current,
							"text",
							lv_text_4_0,
							"org.eclipse.xtext.common.Terminals.STRING");
					}
				)
			)
		)
		    |
		(
			{
				newCompositeNode(grammarAccess.getLabelElementAccess().getKeyLabelsParserRuleCall_3_0());
			}
			ruleKeyLabels
			{
				afterParserOrEnumRuleCall();
			}
			otherlv_6=':'
			{
				newLeafNode(otherlv_6, grammarAccess.getLabelElementAccess().getColonKeyword_3_1());
			}
			{
				if ($current==null) {
					$current = createModelElement(grammarAccess.getLabelElementRule());
				}
				newCompositeNode(grammarAccess.getLabelElementAccess().getElkGraphElementLabelsParserRuleCall_3_2());
			}
			this_ElkGraphElementLabels_7=ruleElkGraphElementLabels[$current]
			{
				$current = $this_ElkGraphElementLabels_7.current;
				afterParserOrEnumRuleCall();
			}
		)
		    |
		(
			{
				newCompositeNode(grammarAccess.getLabelElementAccess().getKeyLayoutOptionsParserRuleCall_4_0());
			}
			ruleKeyLayoutOptions
			{
				afterParserOrEnumRuleCall();
			}
			otherlv_9=':'
			{
				newLeafNode(otherlv_9, grammarAccess.getLabelElementAccess().getColonKeyword_4_1());
			}
			{
				if ($current==null) {
					$current = createModelElement(grammarAccess.getLabelElementRule());
				}
				newCompositeNode(grammarAccess.getLabelElementAccess().getElkGraphElementPropertiesParserRuleCall_4_2());
			}
			this_ElkGraphElementProperties_10=ruleElkGraphElementProperties[$current]
			{
				$current = $this_ElkGraphElementProperties_10.current;
				afterParserOrEnumRuleCall();
			}
		)
		    |
		{
			newCompositeNode(grammarAccess.getLabelElementAccess().getJsonMemberParserRuleCall_5());
		}
		ruleJsonMember
		{
			afterParserOrEnumRuleCall();
		}
	)
;

// Entry rule entryRuleElkEdge
entryRuleElkEdge returns [EObject current=null]:
	{ newCompositeNode(grammarAccess.getElkEdgeRule()); }
	iv_ruleElkEdge=ruleElkEdge
	{ $current=$iv_ruleElkEdge.current; }
	EOF;

// Rule ElkEdge
ruleElkEdge returns [EObject current=null]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		otherlv_0='{'
		{
			newLeafNode(otherlv_0, grammarAccess.getElkEdgeAccess().getLeftCurlyBracketKeyword_0());
		}
		{
			if ($current==null) {
				$current = createModelElement(grammarAccess.getElkEdgeRule());
			}
			newCompositeNode(grammarAccess.getElkEdgeAccess().getEdgeElementParserRuleCall_1());
		}
		this_EdgeElement_1=ruleEdgeElement[$current]
		{
			$current = $this_EdgeElement_1.current;
			afterParserOrEnumRuleCall();
		}
		(
			otherlv_2=','
			{
				newLeafNode(otherlv_2, grammarAccess.getElkEdgeAccess().getCommaKeyword_2_0());
			}
			{
				if ($current==null) {
					$current = createModelElement(grammarAccess.getElkEdgeRule());
				}
				newCompositeNode(grammarAccess.getElkEdgeAccess().getEdgeElementParserRuleCall_2_1());
			}
			this_EdgeElement_3=ruleEdgeElement[$current]
			{
				$current = $this_EdgeElement_3.current;
				afterParserOrEnumRuleCall();
			}
		)*
		(
			otherlv_4=','
			{
				newLeafNode(otherlv_4, grammarAccess.getElkEdgeAccess().getCommaKeyword_3());
			}
		)?
		otherlv_5='}'
		{
			newLeafNode(otherlv_5, grammarAccess.getElkEdgeAccess().getRightCurlyBracketKeyword_4());
		}
	)
;


// Rule EdgeElement
ruleEdgeElement[EObject in_current]  returns [EObject current=in_current]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		{
			if ($current==null) {
				$current = createModelElement(grammarAccess.getEdgeElementRule());
			}
			newCompositeNode(grammarAccess.getEdgeElementAccess().getElkIdParserRuleCall_0());
		}
		this_ElkId_0=ruleElkId[$current]
		{
			$current = $this_ElkId_0.current;
			afterParserOrEnumRuleCall();
		}
		    |
		(
			{
				newCompositeNode(grammarAccess.getEdgeElementAccess().getKeySourcesParserRuleCall_1_0());
			}
			ruleKeySources
			{
				afterParserOrEnumRuleCall();
			}
			otherlv_2=':'
			{
				newLeafNode(otherlv_2, grammarAccess.getEdgeElementAccess().getColonKeyword_1_1());
			}
			{
				if ($current==null) {
					$current = createModelElement(grammarAccess.getEdgeElementRule());
				}
				newCompositeNode(grammarAccess.getEdgeElementAccess().getElkEdgeSourcesParserRuleCall_1_2());
			}
			this_ElkEdgeSources_3=ruleElkEdgeSources[$current]
			{
				$current = $this_ElkEdgeSources_3.current;
				afterParserOrEnumRuleCall();
			}
		)
		    |
		(
			{
				newCompositeNode(grammarAccess.getEdgeElementAccess().getKeyTargetsParserRuleCall_2_0());
			}
			ruleKeyTargets
			{
				afterParserOrEnumRuleCall();
			}
			otherlv_5=':'
			{
				newLeafNode(otherlv_5, grammarAccess.getEdgeElementAccess().getColonKeyword_2_1());
			}
			{
				if ($current==null) {
					$current = createModelElement(grammarAccess.getEdgeElementRule());
				}
				newCompositeNode(grammarAccess.getEdgeElementAccess().getElkEdgeTargetsParserRuleCall_2_2());
			}
			this_ElkEdgeTargets_6=ruleElkEdgeTargets[$current]
			{
				$current = $this_ElkEdgeTargets_6.current;
				afterParserOrEnumRuleCall();
			}
		)
		    |
		(
			{
				newCompositeNode(grammarAccess.getEdgeElementAccess().getKeyLabelsParserRuleCall_3_0());
			}
			ruleKeyLabels
			{
				afterParserOrEnumRuleCall();
			}
			otherlv_8=':'
			{
				newLeafNode(otherlv_8, grammarAccess.getEdgeElementAccess().getColonKeyword_3_1());
			}
			{
				if ($current==null) {
					$current = createModelElement(grammarAccess.getEdgeElementRule());
				}
				newCompositeNode(grammarAccess.getEdgeElementAccess().getElkGraphElementLabelsParserRuleCall_3_2());
			}
			this_ElkGraphElementLabels_9=ruleElkGraphElementLabels[$current]
			{
				$current = $this_ElkGraphElementLabels_9.current;
				afterParserOrEnumRuleCall();
			}
		)
		    |
		(
			{
				newCompositeNode(grammarAccess.getEdgeElementAccess().getKeyLayoutOptionsParserRuleCall_4_0());
			}
			ruleKeyLayoutOptions
			{
				afterParserOrEnumRuleCall();
			}
			otherlv_11=':'
			{
				newLeafNode(otherlv_11, grammarAccess.getEdgeElementAccess().getColonKeyword_4_1());
			}
			{
				if ($current==null) {
					$current = createModelElement(grammarAccess.getEdgeElementRule());
				}
				newCompositeNode(grammarAccess.getEdgeElementAccess().getElkGraphElementPropertiesParserRuleCall_4_2());
			}
			this_ElkGraphElementProperties_12=ruleElkGraphElementProperties[$current]
			{
				$current = $this_ElkGraphElementProperties_12.current;
				afterParserOrEnumRuleCall();
			}
		)
		    |
		{
			newCompositeNode(grammarAccess.getEdgeElementAccess().getJsonMemberParserRuleCall_5());
		}
		ruleJsonMember
		{
			afterParserOrEnumRuleCall();
		}
	)
;


// Rule ElkEdgeSources
ruleElkEdgeSources[EObject in_current]  returns [EObject current=in_current]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		otherlv_0='['
		{
			newLeafNode(otherlv_0, grammarAccess.getElkEdgeSourcesAccess().getLeftSquareBracketKeyword_0());
		}
		(
			(
				(
					{
						if ($current==null) {
							$current = createModelElement(grammarAccess.getElkEdgeSourcesRule());
						}
					}
					otherlv_1=RULE_STRING
					{
						newLeafNode(otherlv_1, grammarAccess.getElkEdgeSourcesAccess().getSourcesElkConnectableShapeCrossReference_1_0_0());
					}
				)
			)
			(
				otherlv_2=','
				{
					newLeafNode(otherlv_2, grammarAccess.getElkEdgeSourcesAccess().getCommaKeyword_1_1_0());
				}
				(
					(
						{
							if ($current==null) {
								$current = createModelElement(grammarAccess.getElkEdgeSourcesRule());
							}
						}
						otherlv_3=RULE_STRING
						{
							newLeafNode(otherlv_3, grammarAccess.getElkEdgeSourcesAccess().getSourcesElkConnectableShapeCrossReference_1_1_1_0());
						}
					)
				)
			)*
		)?
		(
			otherlv_4=','
			{
				newLeafNode(otherlv_4, grammarAccess.getElkEdgeSourcesAccess().getCommaKeyword_2());
			}
		)?
		otherlv_5=']'
		{
			newLeafNode(otherlv_5, grammarAccess.getElkEdgeSourcesAccess().getRightSquareBracketKeyword_3());
		}
	)
;


// Rule ElkEdgeTargets
ruleElkEdgeTargets[EObject in_current]  returns [EObject current=in_current]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		otherlv_0='['
		{
			newLeafNode(otherlv_0, grammarAccess.getElkEdgeTargetsAccess().getLeftSquareBracketKeyword_0());
		}
		(
			(
				(
					{
						if ($current==null) {
							$current = createModelElement(grammarAccess.getElkEdgeTargetsRule());
						}
					}
					otherlv_1=RULE_STRING
					{
						newLeafNode(otherlv_1, grammarAccess.getElkEdgeTargetsAccess().getTargetsElkConnectableShapeCrossReference_1_0_0());
					}
				)
			)
			(
				otherlv_2=','
				{
					newLeafNode(otherlv_2, grammarAccess.getElkEdgeTargetsAccess().getCommaKeyword_1_1_0());
				}
				(
					(
						{
							if ($current==null) {
								$current = createModelElement(grammarAccess.getElkEdgeTargetsRule());
							}
						}
						otherlv_3=RULE_STRING
						{
							newLeafNode(otherlv_3, grammarAccess.getElkEdgeTargetsAccess().getTargetsElkConnectableShapeCrossReference_1_1_1_0());
						}
					)
				)
			)*
		)?
		(
			otherlv_4=','
			{
				newLeafNode(otherlv_4, grammarAccess.getElkEdgeTargetsAccess().getCommaKeyword_2());
			}
		)?
		otherlv_5=']'
		{
			newLeafNode(otherlv_5, grammarAccess.getElkEdgeTargetsAccess().getRightSquareBracketKeyword_3());
		}
	)
;


// Rule ElkId
ruleElkId[EObject in_current]  returns [EObject current=in_current]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		{
			newCompositeNode(grammarAccess.getElkIdAccess().getKeyIdParserRuleCall_0());
		}
		ruleKeyId
		{
			afterParserOrEnumRuleCall();
		}
		otherlv_1=':'
		{
			newLeafNode(otherlv_1, grammarAccess.getElkIdAccess().getColonKeyword_1());
		}
		(
			(
				lv_identifier_2_0=RULE_STRING
				{
					newLeafNode(lv_identifier_2_0, grammarAccess.getElkIdAccess().getIdentifierSTRINGTerminalRuleCall_2_0());
				}
				{
					if ($current==null) {
						$current = createModelElement(grammarAccess.getElkIdRule());
					}
					setWithLastConsumed(
						$current,
						"identifier",
						lv_identifier_2_0,
						"org.eclipse.xtext.common.Terminals.STRING");
				}
			)
		)
	)
;


// Rule ElkNodeChildren
ruleElkNodeChildren[EObject in_current]  returns [EObject current=in_current]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		otherlv_0='['
		{
			newLeafNode(otherlv_0, grammarAccess.getElkNodeChildrenAccess().getLeftSquareBracketKeyword_0());
		}
		(
			(
				(
					{
						newCompositeNode(grammarAccess.getElkNodeChildrenAccess().getChildrenElkNodeParserRuleCall_1_0_0());
					}
					lv_children_1_0=ruleElkNode
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getElkNodeChildrenRule());
						}
						add(
							$current,
							"children",
							lv_children_1_0,
							"org.eclipse.elk.graph.json.text.ElkGraphJson.ElkNode");
						afterParserOrEnumRuleCall();
					}
				)
			)
			(
				otherlv_2=','
				{
					newLeafNode(otherlv_2, grammarAccess.getElkNodeChildrenAccess().getCommaKeyword_1_1_0());
				}
				(
					(
						{
							newCompositeNode(grammarAccess.getElkNodeChildrenAccess().getChildrenElkNodeParserRuleCall_1_1_1_0());
						}
						lv_children_3_0=ruleElkNode
						{
							if ($current==null) {
								$current = createModelElementForParent(grammarAccess.getElkNodeChildrenRule());
							}
							add(
								$current,
								"children",
								lv_children_3_0,
								"org.eclipse.elk.graph.json.text.ElkGraphJson.ElkNode");
							afterParserOrEnumRuleCall();
						}
					)
				)
			)*
		)?
		(
			otherlv_4=','
			{
				newLeafNode(otherlv_4, grammarAccess.getElkNodeChildrenAccess().getCommaKeyword_2());
			}
		)?
		otherlv_5=']'
		{
			newLeafNode(otherlv_5, grammarAccess.getElkNodeChildrenAccess().getRightSquareBracketKeyword_3());
		}
	)
;


// Rule ElkNodePorts
ruleElkNodePorts[EObject in_current]  returns [EObject current=in_current]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		otherlv_0='['
		{
			newLeafNode(otherlv_0, grammarAccess.getElkNodePortsAccess().getLeftSquareBracketKeyword_0());
		}
		(
			(
				(
					{
						newCompositeNode(grammarAccess.getElkNodePortsAccess().getPortsElkPortParserRuleCall_1_0_0());
					}
					lv_ports_1_0=ruleElkPort
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getElkNodePortsRule());
						}
						add(
							$current,
							"ports",
							lv_ports_1_0,
							"org.eclipse.elk.graph.json.text.ElkGraphJson.ElkPort");
						afterParserOrEnumRuleCall();
					}
				)
			)
			(
				otherlv_2=','
				{
					newLeafNode(otherlv_2, grammarAccess.getElkNodePortsAccess().getCommaKeyword_1_1_0());
				}
				(
					(
						{
							newCompositeNode(grammarAccess.getElkNodePortsAccess().getPortsElkPortParserRuleCall_1_1_1_0());
						}
						lv_ports_3_0=ruleElkPort
						{
							if ($current==null) {
								$current = createModelElementForParent(grammarAccess.getElkNodePortsRule());
							}
							add(
								$current,
								"ports",
								lv_ports_3_0,
								"org.eclipse.elk.graph.json.text.ElkGraphJson.ElkPort");
							afterParserOrEnumRuleCall();
						}
					)
				)
			)*
		)?
		(
			otherlv_4=','
			{
				newLeafNode(otherlv_4, grammarAccess.getElkNodePortsAccess().getCommaKeyword_2());
			}
		)?
		otherlv_5=']'
		{
			newLeafNode(otherlv_5, grammarAccess.getElkNodePortsAccess().getRightSquareBracketKeyword_3());
		}
	)
;


// Rule ElkNodeEdges
ruleElkNodeEdges[EObject in_current]  returns [EObject current=in_current]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		otherlv_0='['
		{
			newLeafNode(otherlv_0, grammarAccess.getElkNodeEdgesAccess().getLeftSquareBracketKeyword_0());
		}
		(
			(
				(
					{
						newCompositeNode(grammarAccess.getElkNodeEdgesAccess().getContainedEdgesElkEdgeParserRuleCall_1_0_0());
					}
					lv_containedEdges_1_0=ruleElkEdge
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getElkNodeEdgesRule());
						}
						add(
							$current,
							"containedEdges",
							lv_containedEdges_1_0,
							"org.eclipse.elk.graph.json.text.ElkGraphJson.ElkEdge");
						afterParserOrEnumRuleCall();
					}
				)
			)
			(
				otherlv_2=','
				{
					newLeafNode(otherlv_2, grammarAccess.getElkNodeEdgesAccess().getCommaKeyword_1_1_0());
				}
				(
					(
						{
							newCompositeNode(grammarAccess.getElkNodeEdgesAccess().getContainedEdgesElkEdgeParserRuleCall_1_1_1_0());
						}
						lv_containedEdges_3_0=ruleElkEdge
						{
							if ($current==null) {
								$current = createModelElementForParent(grammarAccess.getElkNodeEdgesRule());
							}
							add(
								$current,
								"containedEdges",
								lv_containedEdges_3_0,
								"org.eclipse.elk.graph.json.text.ElkGraphJson.ElkEdge");
							afterParserOrEnumRuleCall();
						}
					)
				)
			)*
		)?
		(
			otherlv_4=','
			{
				newLeafNode(otherlv_4, grammarAccess.getElkNodeEdgesAccess().getCommaKeyword_2());
			}
		)?
		otherlv_5=']'
		{
			newLeafNode(otherlv_5, grammarAccess.getElkNodeEdgesAccess().getRightSquareBracketKeyword_3());
		}
	)
;


// Rule ElkGraphElementLabels
ruleElkGraphElementLabels[EObject in_current]  returns [EObject current=in_current]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		otherlv_0='['
		{
			newLeafNode(otherlv_0, grammarAccess.getElkGraphElementLabelsAccess().getLeftSquareBracketKeyword_0());
		}
		(
			(
				(
					{
						newCompositeNode(grammarAccess.getElkGraphElementLabelsAccess().getLabelsElkLabelParserRuleCall_1_0_0());
					}
					lv_labels_1_0=ruleElkLabel
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getElkGraphElementLabelsRule());
						}
						add(
							$current,
							"labels",
							lv_labels_1_0,
							"org.eclipse.elk.graph.json.text.ElkGraphJson.ElkLabel");
						afterParserOrEnumRuleCall();
					}
				)
			)
			(
				otherlv_2=','
				{
					newLeafNode(otherlv_2, grammarAccess.getElkGraphElementLabelsAccess().getCommaKeyword_1_1_0());
				}
				(
					(
						{
							newCompositeNode(grammarAccess.getElkGraphElementLabelsAccess().getLabelsElkLabelParserRuleCall_1_1_1_0());
						}
						lv_labels_3_0=ruleElkLabel
						{
							if ($current==null) {
								$current = createModelElementForParent(grammarAccess.getElkGraphElementLabelsRule());
							}
							add(
								$current,
								"labels",
								lv_labels_3_0,
								"org.eclipse.elk.graph.json.text.ElkGraphJson.ElkLabel");
							afterParserOrEnumRuleCall();
						}
					)
				)
			)*
		)?
		(
			otherlv_4=','
			{
				newLeafNode(otherlv_4, grammarAccess.getElkGraphElementLabelsAccess().getCommaKeyword_2());
			}
		)?
		otherlv_5=']'
		{
			newLeafNode(otherlv_5, grammarAccess.getElkGraphElementLabelsAccess().getRightSquareBracketKeyword_3());
		}
	)
;


// Rule ElkGraphElementProperties
ruleElkGraphElementProperties[EObject in_current]  returns [EObject current=in_current]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		otherlv_0='{'
		{
			newLeafNode(otherlv_0, grammarAccess.getElkGraphElementPropertiesAccess().getLeftCurlyBracketKeyword_0());
		}
		(
			(
				(
					{
						newCompositeNode(grammarAccess.getElkGraphElementPropertiesAccess().getPropertiesPropertyParserRuleCall_1_0_0());
					}
					lv_properties_1_0=ruleProperty
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getElkGraphElementPropertiesRule());
						}
						add(
							$current,
							"properties",
							lv_properties_1_0,
							"org.eclipse.elk.graph.json.text.ElkGraphJson.Property");
						afterParserOrEnumRuleCall();
					}
				)
			)
			(
				otherlv_2=','
				{
					newLeafNode(otherlv_2, grammarAccess.getElkGraphElementPropertiesAccess().getCommaKeyword_1_1_0());
				}
				(
					(
						{
							newCompositeNode(grammarAccess.getElkGraphElementPropertiesAccess().getPropertiesPropertyParserRuleCall_1_1_1_0());
						}
						lv_properties_3_0=ruleProperty
						{
							if ($current==null) {
								$current = createModelElementForParent(grammarAccess.getElkGraphElementPropertiesRule());
							}
							add(
								$current,
								"properties",
								lv_properties_3_0,
								"org.eclipse.elk.graph.json.text.ElkGraphJson.Property");
							afterParserOrEnumRuleCall();
						}
					)
				)
			)*
		)?
		(
			otherlv_4=','
			{
				newLeafNode(otherlv_4, grammarAccess.getElkGraphElementPropertiesAccess().getCommaKeyword_2());
			}
		)?
		otherlv_5='}'
		{
			newLeafNode(otherlv_5, grammarAccess.getElkGraphElementPropertiesAccess().getRightCurlyBracketKeyword_3());
		}
	)
;


// Rule ShapeElement
ruleShapeElement[EObject in_current]  returns [EObject current=in_current]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		(
			{
				newCompositeNode(grammarAccess.getShapeElementAccess().getKeyXParserRuleCall_0_0());
			}
			ruleKeyX
			{
				afterParserOrEnumRuleCall();
			}
			otherlv_1=':'
			{
				newLeafNode(otherlv_1, grammarAccess.getShapeElementAccess().getColonKeyword_0_1());
			}
			(
				(
					{
						newCompositeNode(grammarAccess.getShapeElementAccess().getXNumberParserRuleCall_0_2_0());
					}
					lv_x_2_0=ruleNumber
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getShapeElementRule());
						}
						set(
							$current,
							"x",
							lv_x_2_0,
							"org.eclipse.elk.graph.json.text.ElkGraphJson.Number");
						afterParserOrEnumRuleCall();
					}
				)
			)
		)
		    |
		(
			{
				newCompositeNode(grammarAccess.getShapeElementAccess().getKeyYParserRuleCall_1_0());
			}
			ruleKeyY
			{
				afterParserOrEnumRuleCall();
			}
			otherlv_4=':'
			{
				newLeafNode(otherlv_4, grammarAccess.getShapeElementAccess().getColonKeyword_1_1());
			}
			(
				(
					{
						newCompositeNode(grammarAccess.getShapeElementAccess().getYNumberParserRuleCall_1_2_0());
					}
					lv_y_5_0=ruleNumber
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getShapeElementRule());
						}
						set(
							$current,
							"y",
							lv_y_5_0,
							"org.eclipse.elk.graph.json.text.ElkGraphJson.Number");
						afterParserOrEnumRuleCall();
					}
				)
			)
		)
		    |
		(
			{
				newCompositeNode(grammarAccess.getShapeElementAccess().getKeyWidthParserRuleCall_2_0());
			}
			ruleKeyWidth
			{
				afterParserOrEnumRuleCall();
			}
			otherlv_7=':'
			{
				newLeafNode(otherlv_7, grammarAccess.getShapeElementAccess().getColonKeyword_2_1());
			}
			(
				(
					{
						newCompositeNode(grammarAccess.getShapeElementAccess().getWidthNumberParserRuleCall_2_2_0());
					}
					lv_width_8_0=ruleNumber
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getShapeElementRule());
						}
						set(
							$current,
							"width",
							lv_width_8_0,
							"org.eclipse.elk.graph.json.text.ElkGraphJson.Number");
						afterParserOrEnumRuleCall();
					}
				)
			)
		)
		    |
		(
			{
				newCompositeNode(grammarAccess.getShapeElementAccess().getKeyHeightParserRuleCall_3_0());
			}
			ruleKeyHeight
			{
				afterParserOrEnumRuleCall();
			}
			otherlv_10=':'
			{
				newLeafNode(otherlv_10, grammarAccess.getShapeElementAccess().getColonKeyword_3_1());
			}
			(
				(
					{
						newCompositeNode(grammarAccess.getShapeElementAccess().getHeightNumberParserRuleCall_3_2_0());
					}
					lv_height_11_0=ruleNumber
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getShapeElementRule());
						}
						set(
							$current,
							"height",
							lv_height_11_0,
							"org.eclipse.elk.graph.json.text.ElkGraphJson.Number");
						afterParserOrEnumRuleCall();
					}
				)
			)
		)
	)
;

// Entry rule entryRuleProperty
entryRuleProperty returns [EObject current=null]:
	{ newCompositeNode(grammarAccess.getPropertyRule()); }
	iv_ruleProperty=ruleProperty
	{ $current=$iv_ruleProperty.current; }
	EOF;

// Rule Property
ruleProperty returns [EObject current=null]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		(
			(
				{
					newCompositeNode(grammarAccess.getPropertyAccess().getKeyPropertyKeyParserRuleCall_0_0());
				}
				lv_key_0_0=rulePropertyKey
				{
					if ($current==null) {
						$current = createModelElementForParent(grammarAccess.getPropertyRule());
					}
					set(
						$current,
						"key",
						lv_key_0_0,
						"org.eclipse.elk.graph.json.text.ElkGraphJson.PropertyKey");
					afterParserOrEnumRuleCall();
				}
			)
		)
		otherlv_1=':'
		{
			newLeafNode(otherlv_1, grammarAccess.getPropertyAccess().getColonKeyword_1());
		}
		(
			(
				(
					{
						newCompositeNode(grammarAccess.getPropertyAccess().getValueStringValueParserRuleCall_2_0_0());
					}
					lv_value_2_0=ruleStringValue
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getPropertyRule());
						}
						set(
							$current,
							"value",
							lv_value_2_0,
							"org.eclipse.elk.graph.json.text.ElkGraphJson.StringValue");
						afterParserOrEnumRuleCall();
					}
				)
			)
			    |
			(
				(
					{
						newCompositeNode(grammarAccess.getPropertyAccess().getValueNumberValueParserRuleCall_2_1_0());
					}
					lv_value_3_0=ruleNumberValue
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getPropertyRule());
						}
						set(
							$current,
							"value",
							lv_value_3_0,
							"org.eclipse.elk.graph.json.text.ElkGraphJson.NumberValue");
						afterParserOrEnumRuleCall();
					}
				)
			)
			    |
			(
				(
					{
						newCompositeNode(grammarAccess.getPropertyAccess().getValueBooleanValueParserRuleCall_2_2_0());
					}
					lv_value_4_0=ruleBooleanValue
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getPropertyRule());
						}
						set(
							$current,
							"value",
							lv_value_4_0,
							"org.eclipse.elk.graph.json.text.ElkGraphJson.BooleanValue");
						afterParserOrEnumRuleCall();
					}
				)
			)
		)
	)
;

// Entry rule entryRulePropertyKey
entryRulePropertyKey returns [String current=null]@init {
	HiddenTokens myHiddenTokenState = ((XtextTokenStream)input).setHiddenTokens();
}:
	{ newCompositeNode(grammarAccess.getPropertyKeyRule()); }
	iv_rulePropertyKey=rulePropertyKey
	{ $current=$iv_rulePropertyKey.current.getText(); }
	EOF;
finally {
	myHiddenTokenState.restore();
}

// Rule PropertyKey
rulePropertyKey returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()]
@init {
	enterRule();
	HiddenTokens myHiddenTokenState = ((XtextTokenStream)input).setHiddenTokens();
}
@after {
	leaveRule();
}:
	(
		this_STRING_0=RULE_STRING
		{
			$current.merge(this_STRING_0);
		}
		{
			newLeafNode(this_STRING_0, grammarAccess.getPropertyKeyAccess().getSTRINGTerminalRuleCall_0());
		}
		    |
		this_ID_1=RULE_ID
		{
			$current.merge(this_ID_1);
		}
		{
			newLeafNode(this_ID_1, grammarAccess.getPropertyKeyAccess().getIDTerminalRuleCall_1());
		}
	)
;
finally {
	myHiddenTokenState.restore();
}

// Entry rule entryRuleStringValue
entryRuleStringValue returns [String current=null]:
	{ newCompositeNode(grammarAccess.getStringValueRule()); }
	iv_ruleStringValue=ruleStringValue
	{ $current=$iv_ruleStringValue.current.getText(); }
	EOF;

// Rule StringValue
ruleStringValue returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	this_STRING_0=RULE_STRING
	{
		$current.merge(this_STRING_0);
	}
	{
		newLeafNode(this_STRING_0, grammarAccess.getStringValueAccess().getSTRINGTerminalRuleCall());
	}
;

// Entry rule entryRuleNumberValue
entryRuleNumberValue returns [String current=null]:
	{ newCompositeNode(grammarAccess.getNumberValueRule()); }
	iv_ruleNumberValue=ruleNumberValue
	{ $current=$iv_ruleNumberValue.current.getText(); }
	EOF;

// Rule NumberValue
ruleNumberValue returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		this_SIGNED_INT_0=RULE_SIGNED_INT
		{
			$current.merge(this_SIGNED_INT_0);
		}
		{
			newLeafNode(this_SIGNED_INT_0, grammarAccess.getNumberValueAccess().getSIGNED_INTTerminalRuleCall_0());
		}
		    |
		this_FLOAT_1=RULE_FLOAT
		{
			$current.merge(this_FLOAT_1);
		}
		{
			newLeafNode(this_FLOAT_1, grammarAccess.getNumberValueAccess().getFLOATTerminalRuleCall_1());
		}
	)
;

// Entry rule entryRuleBooleanValue
entryRuleBooleanValue returns [String current=null]:
	{ newCompositeNode(grammarAccess.getBooleanValueRule()); }
	iv_ruleBooleanValue=ruleBooleanValue
	{ $current=$iv_ruleBooleanValue.current.getText(); }
	EOF;

// Rule BooleanValue
ruleBooleanValue returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		kw='true'
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getBooleanValueAccess().getTrueKeyword_0());
		}
		    |
		kw='false'
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getBooleanValueAccess().getFalseKeyword_1());
		}
	)
;

// Entry rule entryRuleNumber
entryRuleNumber returns [String current=null]:
	{ newCompositeNode(grammarAccess.getNumberRule()); }
	iv_ruleNumber=ruleNumber
	{ $current=$iv_ruleNumber.current.getText(); }
	EOF;

// Rule Number
ruleNumber returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		this_SIGNED_INT_0=RULE_SIGNED_INT
		{
			$current.merge(this_SIGNED_INT_0);
		}
		{
			newLeafNode(this_SIGNED_INT_0, grammarAccess.getNumberAccess().getSIGNED_INTTerminalRuleCall_0());
		}
		    |
		this_FLOAT_1=RULE_FLOAT
		{
			$current.merge(this_FLOAT_1);
		}
		{
			newLeafNode(this_FLOAT_1, grammarAccess.getNumberAccess().getFLOATTerminalRuleCall_1());
		}
	)
;

// Entry rule entryRuleJsonObject
entryRuleJsonObject returns [String current=null]:
	{ newCompositeNode(grammarAccess.getJsonObjectRule()); }
	iv_ruleJsonObject=ruleJsonObject
	{ $current=$iv_ruleJsonObject.current.getText(); }
	EOF;

// Rule JsonObject
ruleJsonObject returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		kw='{'
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getJsonObjectAccess().getLeftCurlyBracketKeyword_0());
		}
		(
			{
				newCompositeNode(grammarAccess.getJsonObjectAccess().getJsonMemberParserRuleCall_1_0());
			}
			this_JsonMember_1=ruleJsonMember
			{
				$current.merge(this_JsonMember_1);
			}
			{
				afterParserOrEnumRuleCall();
			}
			(
				kw=','
				{
					$current.merge(kw);
					newLeafNode(kw, grammarAccess.getJsonObjectAccess().getCommaKeyword_1_1_0());
				}
				{
					newCompositeNode(grammarAccess.getJsonObjectAccess().getJsonMemberParserRuleCall_1_1_1());
				}
				this_JsonMember_3=ruleJsonMember
				{
					$current.merge(this_JsonMember_3);
				}
				{
					afterParserOrEnumRuleCall();
				}
			)*
		)?
		(
			kw=','
			{
				$current.merge(kw);
				newLeafNode(kw, grammarAccess.getJsonObjectAccess().getCommaKeyword_2());
			}
		)?
		kw='}'
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getJsonObjectAccess().getRightCurlyBracketKeyword_3());
		}
	)
;

// Entry rule entryRuleJsonArray
entryRuleJsonArray returns [String current=null]:
	{ newCompositeNode(grammarAccess.getJsonArrayRule()); }
	iv_ruleJsonArray=ruleJsonArray
	{ $current=$iv_ruleJsonArray.current.getText(); }
	EOF;

// Rule JsonArray
ruleJsonArray returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		kw='['
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getJsonArrayAccess().getLeftSquareBracketKeyword_0());
		}
		(
			{
				newCompositeNode(grammarAccess.getJsonArrayAccess().getJsonValueParserRuleCall_1_0());
			}
			this_JsonValue_1=ruleJsonValue
			{
				$current.merge(this_JsonValue_1);
			}
			{
				afterParserOrEnumRuleCall();
			}
			(
				kw=','
				{
					$current.merge(kw);
					newLeafNode(kw, grammarAccess.getJsonArrayAccess().getCommaKeyword_1_1_0());
				}
				{
					newCompositeNode(grammarAccess.getJsonArrayAccess().getJsonValueParserRuleCall_1_1_1());
				}
				this_JsonValue_3=ruleJsonValue
				{
					$current.merge(this_JsonValue_3);
				}
				{
					afterParserOrEnumRuleCall();
				}
			)*
		)?
		(
			kw=','
			{
				$current.merge(kw);
				newLeafNode(kw, grammarAccess.getJsonArrayAccess().getCommaKeyword_2());
			}
		)?
		kw=']'
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getJsonArrayAccess().getRightSquareBracketKeyword_3());
		}
	)
;

// Entry rule entryRuleJsonMember
entryRuleJsonMember returns [String current=null]:
	{ newCompositeNode(grammarAccess.getJsonMemberRule()); }
	iv_ruleJsonMember=ruleJsonMember
	{ $current=$iv_ruleJsonMember.current.getText(); }
	EOF;

// Rule JsonMember
ruleJsonMember returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		(
			this_STRING_0=RULE_STRING
			{
				$current.merge(this_STRING_0);
			}
			{
				newLeafNode(this_STRING_0, grammarAccess.getJsonMemberAccess().getSTRINGTerminalRuleCall_0_0());
			}
			    |
			this_ID_1=RULE_ID
			{
				$current.merge(this_ID_1);
			}
			{
				newLeafNode(this_ID_1, grammarAccess.getJsonMemberAccess().getIDTerminalRuleCall_0_1());
			}
		)
		kw=':'
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getJsonMemberAccess().getColonKeyword_1());
		}
		{
			newCompositeNode(grammarAccess.getJsonMemberAccess().getJsonValueParserRuleCall_2());
		}
		this_JsonValue_3=ruleJsonValue
		{
			$current.merge(this_JsonValue_3);
		}
		{
			afterParserOrEnumRuleCall();
		}
	)
;

// Entry rule entryRuleJsonValue
entryRuleJsonValue returns [String current=null]:
	{ newCompositeNode(grammarAccess.getJsonValueRule()); }
	iv_ruleJsonValue=ruleJsonValue
	{ $current=$iv_ruleJsonValue.current.getText(); }
	EOF;

// Rule JsonValue
ruleJsonValue returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		this_STRING_0=RULE_STRING
		{
			$current.merge(this_STRING_0);
		}
		{
			newLeafNode(this_STRING_0, grammarAccess.getJsonValueAccess().getSTRINGTerminalRuleCall_0());
		}
		    |
		{
			newCompositeNode(grammarAccess.getJsonValueAccess().getNumberParserRuleCall_1());
		}
		this_Number_1=ruleNumber
		{
			$current.merge(this_Number_1);
		}
		{
			afterParserOrEnumRuleCall();
		}
		    |
		{
			newCompositeNode(grammarAccess.getJsonValueAccess().getJsonObjectParserRuleCall_2());
		}
		this_JsonObject_2=ruleJsonObject
		{
			$current.merge(this_JsonObject_2);
		}
		{
			afterParserOrEnumRuleCall();
		}
		    |
		{
			newCompositeNode(grammarAccess.getJsonValueAccess().getJsonArrayParserRuleCall_3());
		}
		this_JsonArray_3=ruleJsonArray
		{
			$current.merge(this_JsonArray_3);
		}
		{
			afterParserOrEnumRuleCall();
		}
		    |
		{
			newCompositeNode(grammarAccess.getJsonValueAccess().getBooleanValueParserRuleCall_4());
		}
		this_BooleanValue_4=ruleBooleanValue
		{
			$current.merge(this_BooleanValue_4);
		}
		{
			afterParserOrEnumRuleCall();
		}
		    |
		kw='null'
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getJsonValueAccess().getNullKeyword_5());
		}
	)
;

// Entry rule entryRuleKeyChildren
entryRuleKeyChildren returns [String current=null]:
	{ newCompositeNode(grammarAccess.getKeyChildrenRule()); }
	iv_ruleKeyChildren=ruleKeyChildren
	{ $current=$iv_ruleKeyChildren.current.getText(); }
	EOF;

// Rule KeyChildren
ruleKeyChildren returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		kw='"children"'
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyChildrenAccess().getChildrenKeyword_0());
		}
		    |
		kw='\'children\''
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyChildrenAccess().getChildrenKeyword_1());
		}
		    |
		kw='children'
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyChildrenAccess().getChildrenKeyword_2());
		}
	)
;

// Entry rule entryRuleKeyPorts
entryRuleKeyPorts returns [String current=null]:
	{ newCompositeNode(grammarAccess.getKeyPortsRule()); }
	iv_ruleKeyPorts=ruleKeyPorts
	{ $current=$iv_ruleKeyPorts.current.getText(); }
	EOF;

// Rule KeyPorts
ruleKeyPorts returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		kw='"ports"'
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyPortsAccess().getPortsKeyword_0());
		}
		    |
		kw='\'ports\''
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyPortsAccess().getPortsKeyword_1());
		}
		    |
		kw='ports'
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyPortsAccess().getPortsKeyword_2());
		}
	)
;

// Entry rule entryRuleKeyLabels
entryRuleKeyLabels returns [String current=null]:
	{ newCompositeNode(grammarAccess.getKeyLabelsRule()); }
	iv_ruleKeyLabels=ruleKeyLabels
	{ $current=$iv_ruleKeyLabels.current.getText(); }
	EOF;

// Rule KeyLabels
ruleKeyLabels returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		kw='"labels"'
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyLabelsAccess().getLabelsKeyword_0());
		}
		    |
		kw='\'labels\''
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyLabelsAccess().getLabelsKeyword_1());
		}
		    |
		kw='labels'
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyLabelsAccess().getLabelsKeyword_2());
		}
	)
;

// Entry rule entryRuleKeyEdges
entryRuleKeyEdges returns [String current=null]:
	{ newCompositeNode(grammarAccess.getKeyEdgesRule()); }
	iv_ruleKeyEdges=ruleKeyEdges
	{ $current=$iv_ruleKeyEdges.current.getText(); }
	EOF;

// Rule KeyEdges
ruleKeyEdges returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		kw='"edges"'
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyEdgesAccess().getEdgesKeyword_0());
		}
		    |
		kw='\'edges\''
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyEdgesAccess().getEdgesKeyword_1());
		}
		    |
		kw='edges'
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyEdgesAccess().getEdgesKeyword_2());
		}
	)
;

// Entry rule entryRuleKeyLayoutOptions
entryRuleKeyLayoutOptions returns [String current=null]:
	{ newCompositeNode(grammarAccess.getKeyLayoutOptionsRule()); }
	iv_ruleKeyLayoutOptions=ruleKeyLayoutOptions
	{ $current=$iv_ruleKeyLayoutOptions.current.getText(); }
	EOF;

// Rule KeyLayoutOptions
ruleKeyLayoutOptions returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		kw='"layoutOptions"'
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyLayoutOptionsAccess().getLayoutOptionsKeyword_0());
		}
		    |
		kw='\'layoutOptions\''
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyLayoutOptionsAccess().getLayoutOptionsKeyword_1());
		}
		    |
		kw='layoutOptions'
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyLayoutOptionsAccess().getLayoutOptionsKeyword_2());
		}
		    |
		kw='"properties"'
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyLayoutOptionsAccess().getPropertiesKeyword_3());
		}
		    |
		kw='\'properties\''
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyLayoutOptionsAccess().getPropertiesKeyword_4());
		}
		    |
		kw='properties'
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyLayoutOptionsAccess().getPropertiesKeyword_5());
		}
	)
;

// Entry rule entryRuleKeyId
entryRuleKeyId returns [String current=null]:
	{ newCompositeNode(grammarAccess.getKeyIdRule()); }
	iv_ruleKeyId=ruleKeyId
	{ $current=$iv_ruleKeyId.current.getText(); }
	EOF;

// Rule KeyId
ruleKeyId returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		kw='"id"'
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyIdAccess().getIdKeyword_0());
		}
		    |
		kw='\'id\''
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyIdAccess().getIdKeyword_1());
		}
		    |
		kw='id'
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyIdAccess().getIdKeyword_2());
		}
	)
;

// Entry rule entryRuleKeyX
entryRuleKeyX returns [String current=null]:
	{ newCompositeNode(grammarAccess.getKeyXRule()); }
	iv_ruleKeyX=ruleKeyX
	{ $current=$iv_ruleKeyX.current.getText(); }
	EOF;

// Rule KeyX
ruleKeyX returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		kw='"x"'
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyXAccess().getXKeyword_0());
		}
		    |
		kw='\'x\''
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyXAccess().getXKeyword_1());
		}
		    |
		kw='x'
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyXAccess().getXKeyword_2());
		}
	)
;

// Entry rule entryRuleKeyY
entryRuleKeyY returns [String current=null]:
	{ newCompositeNode(grammarAccess.getKeyYRule()); }
	iv_ruleKeyY=ruleKeyY
	{ $current=$iv_ruleKeyY.current.getText(); }
	EOF;

// Rule KeyY
ruleKeyY returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		kw='"y"'
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyYAccess().getYKeyword_0());
		}
		    |
		kw='\'y\''
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyYAccess().getYKeyword_1());
		}
		    |
		kw='y'
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyYAccess().getYKeyword_2());
		}
	)
;

// Entry rule entryRuleKeyWidth
entryRuleKeyWidth returns [String current=null]:
	{ newCompositeNode(grammarAccess.getKeyWidthRule()); }
	iv_ruleKeyWidth=ruleKeyWidth
	{ $current=$iv_ruleKeyWidth.current.getText(); }
	EOF;

// Rule KeyWidth
ruleKeyWidth returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		kw='"width"'
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyWidthAccess().getWidthKeyword_0());
		}
		    |
		kw='\'width\''
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyWidthAccess().getWidthKeyword_1());
		}
		    |
		kw='width'
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyWidthAccess().getWidthKeyword_2());
		}
	)
;

// Entry rule entryRuleKeyHeight
entryRuleKeyHeight returns [String current=null]:
	{ newCompositeNode(grammarAccess.getKeyHeightRule()); }
	iv_ruleKeyHeight=ruleKeyHeight
	{ $current=$iv_ruleKeyHeight.current.getText(); }
	EOF;

// Rule KeyHeight
ruleKeyHeight returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		kw='"height"'
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyHeightAccess().getHeightKeyword_0());
		}
		    |
		kw='\'height\''
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyHeightAccess().getHeightKeyword_1());
		}
		    |
		kw='height'
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyHeightAccess().getHeightKeyword_2());
		}
	)
;

// Entry rule entryRuleKeySources
entryRuleKeySources returns [String current=null]:
	{ newCompositeNode(grammarAccess.getKeySourcesRule()); }
	iv_ruleKeySources=ruleKeySources
	{ $current=$iv_ruleKeySources.current.getText(); }
	EOF;

// Rule KeySources
ruleKeySources returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		kw='"sources"'
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeySourcesAccess().getSourcesKeyword_0());
		}
		    |
		kw='\'sources\''
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeySourcesAccess().getSourcesKeyword_1());
		}
		    |
		kw='sources'
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeySourcesAccess().getSourcesKeyword_2());
		}
	)
;

// Entry rule entryRuleKeyTargets
entryRuleKeyTargets returns [String current=null]:
	{ newCompositeNode(grammarAccess.getKeyTargetsRule()); }
	iv_ruleKeyTargets=ruleKeyTargets
	{ $current=$iv_ruleKeyTargets.current.getText(); }
	EOF;

// Rule KeyTargets
ruleKeyTargets returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		kw='"targets"'
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyTargetsAccess().getTargetsKeyword_0());
		}
		    |
		kw='\'targets\''
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyTargetsAccess().getTargetsKeyword_1());
		}
		    |
		kw='targets'
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyTargetsAccess().getTargetsKeyword_2());
		}
	)
;

// Entry rule entryRuleKeyText
entryRuleKeyText returns [String current=null]:
	{ newCompositeNode(grammarAccess.getKeyTextRule()); }
	iv_ruleKeyText=ruleKeyText
	{ $current=$iv_ruleKeyText.current.getText(); }
	EOF;

// Rule KeyText
ruleKeyText returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		kw='"text"'
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyTextAccess().getTextKeyword_0());
		}
		    |
		kw='\'text\''
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyTextAccess().getTextKeyword_1());
		}
		    |
		kw='text'
		{
			$current.merge(kw);
			newLeafNode(kw, grammarAccess.getKeyTextAccess().getTextKeyword_2());
		}
	)
;

RULE_SIGNED_INT : ('+'|'-')? RULE_INT;

RULE_FLOAT : ('+'|'-')? (RULE_INT '.' RULE_INT|RULE_INT ('.' RULE_INT)? ('e'|'E') ('+'|'-')? RULE_INT);

RULE_ID : '^'? ('a'..'z'|'A'..'Z'|'_') ('a'..'z'|'A'..'Z'|'_'|'0'..'9')*;

fragment RULE_INT : ('0'..'9')+;

RULE_STRING : ('"' ('\\' .|~(('\\'|'"')))* '"'|'\'' ('\\' .|~(('\\'|'\'')))* '\'');

RULE_ML_COMMENT : '/*' ( options {greedy=false;} : . )*'*/';

RULE_SL_COMMENT : '//' ~(('\n'|'\r'))* ('\r'? '\n')?;

RULE_WS : (' '|'\t'|'\r'|'\n')+;

RULE_ANY_OTHER : .;
