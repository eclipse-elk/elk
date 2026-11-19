/*******************************************************************************
 * Copyright (c) 2016 TypeFox GmbH (http://www.typefox.io) and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
grammar InternalElkGraph;

options {
	superClass=AbstractInternalAntlrParser;
}

@lexer::header {
package org.eclipse.elk.graph.text.parser.antlr.internal;

// Hack: Use our own Lexer superclass by means of import. 
// Currently there is no other way to specify the superclass for the lexer.
import org.eclipse.xtext.parser.antlr.Lexer;
}

@parser::header {
package org.eclipse.elk.graph.text.parser.antlr.internal;

import org.eclipse.xtext.*;
import org.eclipse.xtext.parser.*;
import org.eclipse.xtext.parser.impl.*;
import org.eclipse.emf.ecore.util.EcoreUtil;
import org.eclipse.emf.ecore.EObject;
import org.eclipse.xtext.parser.antlr.AbstractInternalAntlrParser;
import org.eclipse.xtext.parser.antlr.XtextTokenStream;
import org.eclipse.xtext.parser.antlr.XtextTokenStream.HiddenTokens;
import org.eclipse.xtext.parser.antlr.AntlrDatatypeRuleToken;
import org.eclipse.elk.graph.text.services.ElkGraphGrammarAccess;

}

@parser::members {

 	private ElkGraphGrammarAccess grammarAccess;

    public InternalElkGraphParser(TokenStream input, ElkGraphGrammarAccess grammarAccess) {
        this(input);
        this.grammarAccess = grammarAccess;
        registerRules(grammarAccess.getGrammar());
    }

    @Override
    protected String getFirstRuleName() {
    	return "RootNode";
   	}

   	@Override
   	protected ElkGraphGrammarAccess getGrammarAccess() {
   		return grammarAccess;
   	}

}

@rulecatch {
    catch (RecognitionException re) {
        recover(input,re);
        appendSkippedTokens();
    }
}

// Entry rule entryRuleRootNode
entryRuleRootNode returns [EObject current=null]:
	{ newCompositeNode(grammarAccess.getRootNodeRule()); }
	iv_ruleRootNode=ruleRootNode
	{ $current=$iv_ruleRootNode.current; }
	EOF;

// Rule RootNode
ruleRootNode returns [EObject current=null]
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
					grammarAccess.getRootNodeAccess().getElkNodeAction_0(),
					$current);
			}
		)
		(
			otherlv_1='graph'
			{
				newLeafNode(otherlv_1, grammarAccess.getRootNodeAccess().getGraphKeyword_1_0());
			}
			(
				(
					lv_identifier_2_0=RULE_ID
					{
						newLeafNode(lv_identifier_2_0, grammarAccess.getRootNodeAccess().getIdentifierIDTerminalRuleCall_1_1_0());
					}
					{
						if ($current==null) {
							$current = createModelElement(grammarAccess.getRootNodeRule());
						}
						setWithLastConsumed(
							$current,
							"identifier",
							lv_identifier_2_0,
							"org.eclipse.xtext.common.Terminals.ID");
					}
				)
			)
		)?
		(
			{
				if ($current==null) {
					$current = createModelElement(grammarAccess.getRootNodeRule());
				}
				newCompositeNode(grammarAccess.getRootNodeAccess().getShapeLayoutParserRuleCall_2());
			}
			this_ShapeLayout_3=ruleShapeLayout[$current]
			{
				$current = $this_ShapeLayout_3.current;
				afterParserOrEnumRuleCall();
			}
		)?
		(
			(
				{
					newCompositeNode(grammarAccess.getRootNodeAccess().getPropertiesPropertyParserRuleCall_3_0());
				}
				lv_properties_4_0=ruleProperty
				{
					if ($current==null) {
						$current = createModelElementForParent(grammarAccess.getRootNodeRule());
					}
					add(
						$current,
						"properties",
						lv_properties_4_0,
						"org.eclipse.elk.graph.text.ElkGraph.Property");
					afterParserOrEnumRuleCall();
				}
			)
		)*
		(
			(
				(
					{
						newCompositeNode(grammarAccess.getRootNodeAccess().getLabelsElkLabelParserRuleCall_4_0_0());
					}
					lv_labels_5_0=ruleElkLabel
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getRootNodeRule());
						}
						add(
							$current,
							"labels",
							lv_labels_5_0,
							"org.eclipse.elk.graph.text.ElkGraph.ElkLabel");
						afterParserOrEnumRuleCall();
					}
				)
			)
			    |
			(
				(
					{
						newCompositeNode(grammarAccess.getRootNodeAccess().getPortsElkPortParserRuleCall_4_1_0());
					}
					lv_ports_6_0=ruleElkPort
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getRootNodeRule());
						}
						add(
							$current,
							"ports",
							lv_ports_6_0,
							"org.eclipse.elk.graph.text.ElkGraph.ElkPort");
						afterParserOrEnumRuleCall();
					}
				)
			)
			    |
			(
				(
					{
						newCompositeNode(grammarAccess.getRootNodeAccess().getChildrenElkNodeParserRuleCall_4_2_0());
					}
					lv_children_7_0=ruleElkNode
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getRootNodeRule());
						}
						add(
							$current,
							"children",
							lv_children_7_0,
							"org.eclipse.elk.graph.text.ElkGraph.ElkNode");
						afterParserOrEnumRuleCall();
					}
				)
			)
			    |
			(
				(
					{
						newCompositeNode(grammarAccess.getRootNodeAccess().getContainedEdgesElkEdgeParserRuleCall_4_3_0());
					}
					lv_containedEdges_8_0=ruleElkEdge
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getRootNodeRule());
						}
						add(
							$current,
							"containedEdges",
							lv_containedEdges_8_0,
							"org.eclipse.elk.graph.text.ElkGraph.ElkEdge");
						afterParserOrEnumRuleCall();
					}
				)
			)
		)*
	)
;

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
		otherlv_0='node'
		{
			newLeafNode(otherlv_0, grammarAccess.getElkNodeAccess().getNodeKeyword_0());
		}
		(
			(
				lv_identifier_1_0=RULE_ID
				{
					newLeafNode(lv_identifier_1_0, grammarAccess.getElkNodeAccess().getIdentifierIDTerminalRuleCall_1_0());
				}
				{
					if ($current==null) {
						$current = createModelElement(grammarAccess.getElkNodeRule());
					}
					setWithLastConsumed(
						$current,
						"identifier",
						lv_identifier_1_0,
						"org.eclipse.xtext.common.Terminals.ID");
				}
			)
		)
		(
			otherlv_2='{'
			{
				newLeafNode(otherlv_2, grammarAccess.getElkNodeAccess().getLeftCurlyBracketKeyword_2_0());
			}
			(
				{
					if ($current==null) {
						$current = createModelElement(grammarAccess.getElkNodeRule());
					}
					newCompositeNode(grammarAccess.getElkNodeAccess().getShapeLayoutParserRuleCall_2_1());
				}
				this_ShapeLayout_3=ruleShapeLayout[$current]
				{
					$current = $this_ShapeLayout_3.current;
					afterParserOrEnumRuleCall();
				}
			)?
			(
				(
					{
						newCompositeNode(grammarAccess.getElkNodeAccess().getPropertiesPropertyParserRuleCall_2_2_0());
					}
					lv_properties_4_0=ruleProperty
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getElkNodeRule());
						}
						add(
							$current,
							"properties",
							lv_properties_4_0,
							"org.eclipse.elk.graph.text.ElkGraph.Property");
						afterParserOrEnumRuleCall();
					}
				)
			)*
			(
				(
					(
						{
							newCompositeNode(grammarAccess.getElkNodeAccess().getLabelsElkLabelParserRuleCall_2_3_0_0());
						}
						lv_labels_5_0=ruleElkLabel
						{
							if ($current==null) {
								$current = createModelElementForParent(grammarAccess.getElkNodeRule());
							}
							add(
								$current,
								"labels",
								lv_labels_5_0,
								"org.eclipse.elk.graph.text.ElkGraph.ElkLabel");
							afterParserOrEnumRuleCall();
						}
					)
				)
				    |
				(
					(
						{
							newCompositeNode(grammarAccess.getElkNodeAccess().getPortsElkPortParserRuleCall_2_3_1_0());
						}
						lv_ports_6_0=ruleElkPort
						{
							if ($current==null) {
								$current = createModelElementForParent(grammarAccess.getElkNodeRule());
							}
							add(
								$current,
								"ports",
								lv_ports_6_0,
								"org.eclipse.elk.graph.text.ElkGraph.ElkPort");
							afterParserOrEnumRuleCall();
						}
					)
				)
				    |
				(
					(
						{
							newCompositeNode(grammarAccess.getElkNodeAccess().getChildrenElkNodeParserRuleCall_2_3_2_0());
						}
						lv_children_7_0=ruleElkNode
						{
							if ($current==null) {
								$current = createModelElementForParent(grammarAccess.getElkNodeRule());
							}
							add(
								$current,
								"children",
								lv_children_7_0,
								"org.eclipse.elk.graph.text.ElkGraph.ElkNode");
							afterParserOrEnumRuleCall();
						}
					)
				)
				    |
				(
					(
						{
							newCompositeNode(grammarAccess.getElkNodeAccess().getContainedEdgesElkEdgeParserRuleCall_2_3_3_0());
						}
						lv_containedEdges_8_0=ruleElkEdge
						{
							if ($current==null) {
								$current = createModelElementForParent(grammarAccess.getElkNodeRule());
							}
							add(
								$current,
								"containedEdges",
								lv_containedEdges_8_0,
								"org.eclipse.elk.graph.text.ElkGraph.ElkEdge");
							afterParserOrEnumRuleCall();
						}
					)
				)
			)*
			otherlv_9='}'
			{
				newLeafNode(otherlv_9, grammarAccess.getElkNodeAccess().getRightCurlyBracketKeyword_2_4());
			}
		)?
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
		otherlv_0='label'
		{
			newLeafNode(otherlv_0, grammarAccess.getElkLabelAccess().getLabelKeyword_0());
		}
		(
			(
				(
					lv_identifier_1_0=RULE_ID
					{
						newLeafNode(lv_identifier_1_0, grammarAccess.getElkLabelAccess().getIdentifierIDTerminalRuleCall_1_0_0());
					}
					{
						if ($current==null) {
							$current = createModelElement(grammarAccess.getElkLabelRule());
						}
						setWithLastConsumed(
							$current,
							"identifier",
							lv_identifier_1_0,
							"org.eclipse.xtext.common.Terminals.ID");
					}
				)
			)
			otherlv_2=':'
			{
				newLeafNode(otherlv_2, grammarAccess.getElkLabelAccess().getColonKeyword_1_1());
			}
		)?
		(
			(
				lv_text_3_0=RULE_STRING
				{
					newLeafNode(lv_text_3_0, grammarAccess.getElkLabelAccess().getTextSTRINGTerminalRuleCall_2_0());
				}
				{
					if ($current==null) {
						$current = createModelElement(grammarAccess.getElkLabelRule());
					}
					setWithLastConsumed(
						$current,
						"text",
						lv_text_3_0,
						"org.eclipse.xtext.common.Terminals.STRING");
				}
			)
		)
		(
			otherlv_4='{'
			{
				newLeafNode(otherlv_4, grammarAccess.getElkLabelAccess().getLeftCurlyBracketKeyword_3_0());
			}
			(
				{
					if ($current==null) {
						$current = createModelElement(grammarAccess.getElkLabelRule());
					}
					newCompositeNode(grammarAccess.getElkLabelAccess().getShapeLayoutParserRuleCall_3_1());
				}
				this_ShapeLayout_5=ruleShapeLayout[$current]
				{
					$current = $this_ShapeLayout_5.current;
					afterParserOrEnumRuleCall();
				}
			)?
			(
				(
					{
						newCompositeNode(grammarAccess.getElkLabelAccess().getPropertiesPropertyParserRuleCall_3_2_0());
					}
					lv_properties_6_0=ruleProperty
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getElkLabelRule());
						}
						add(
							$current,
							"properties",
							lv_properties_6_0,
							"org.eclipse.elk.graph.text.ElkGraph.Property");
						afterParserOrEnumRuleCall();
					}
				)
			)*
			(
				(
					{
						newCompositeNode(grammarAccess.getElkLabelAccess().getLabelsElkLabelParserRuleCall_3_3_0());
					}
					lv_labels_7_0=ruleElkLabel
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getElkLabelRule());
						}
						add(
							$current,
							"labels",
							lv_labels_7_0,
							"org.eclipse.elk.graph.text.ElkGraph.ElkLabel");
						afterParserOrEnumRuleCall();
					}
				)
			)*
			otherlv_8='}'
			{
				newLeafNode(otherlv_8, grammarAccess.getElkLabelAccess().getRightCurlyBracketKeyword_3_4());
			}
		)?
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
		otherlv_0='port'
		{
			newLeafNode(otherlv_0, grammarAccess.getElkPortAccess().getPortKeyword_0());
		}
		(
			(
				lv_identifier_1_0=RULE_ID
				{
					newLeafNode(lv_identifier_1_0, grammarAccess.getElkPortAccess().getIdentifierIDTerminalRuleCall_1_0());
				}
				{
					if ($current==null) {
						$current = createModelElement(grammarAccess.getElkPortRule());
					}
					setWithLastConsumed(
						$current,
						"identifier",
						lv_identifier_1_0,
						"org.eclipse.xtext.common.Terminals.ID");
				}
			)
		)
		(
			otherlv_2='{'
			{
				newLeafNode(otherlv_2, grammarAccess.getElkPortAccess().getLeftCurlyBracketKeyword_2_0());
			}
			(
				{
					if ($current==null) {
						$current = createModelElement(grammarAccess.getElkPortRule());
					}
					newCompositeNode(grammarAccess.getElkPortAccess().getShapeLayoutParserRuleCall_2_1());
				}
				this_ShapeLayout_3=ruleShapeLayout[$current]
				{
					$current = $this_ShapeLayout_3.current;
					afterParserOrEnumRuleCall();
				}
			)?
			(
				(
					{
						newCompositeNode(grammarAccess.getElkPortAccess().getPropertiesPropertyParserRuleCall_2_2_0());
					}
					lv_properties_4_0=ruleProperty
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getElkPortRule());
						}
						add(
							$current,
							"properties",
							lv_properties_4_0,
							"org.eclipse.elk.graph.text.ElkGraph.Property");
						afterParserOrEnumRuleCall();
					}
				)
			)*
			(
				(
					{
						newCompositeNode(grammarAccess.getElkPortAccess().getLabelsElkLabelParserRuleCall_2_3_0());
					}
					lv_labels_5_0=ruleElkLabel
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getElkPortRule());
						}
						add(
							$current,
							"labels",
							lv_labels_5_0,
							"org.eclipse.elk.graph.text.ElkGraph.ElkLabel");
						afterParserOrEnumRuleCall();
					}
				)
			)*
			otherlv_6='}'
			{
				newLeafNode(otherlv_6, grammarAccess.getElkPortAccess().getRightCurlyBracketKeyword_2_4());
			}
		)?
	)
;


// Rule ShapeLayout
ruleShapeLayout[EObject in_current]  returns [EObject current=in_current]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		otherlv_0='layout'
		{
			newLeafNode(otherlv_0, grammarAccess.getShapeLayoutAccess().getLayoutKeyword_0());
		}
		otherlv_1='['
		{
			newLeafNode(otherlv_1, grammarAccess.getShapeLayoutAccess().getLeftSquareBracketKeyword_1());
		}
		(
			(
				{ 
				  getUnorderedGroupHelper().enter(grammarAccess.getShapeLayoutAccess().getUnorderedGroup_2());
				}
				(
					(
			(
				{getUnorderedGroupHelper().canSelect(grammarAccess.getShapeLayoutAccess().getUnorderedGroup_2(), 0)}?=>(
					{
						getUnorderedGroupHelper().select(grammarAccess.getShapeLayoutAccess().getUnorderedGroup_2(), 0);
					}
								({true}?=>(otherlv_3='position'
								{
									newLeafNode(otherlv_3, grammarAccess.getShapeLayoutAccess().getPositionKeyword_2_0_0());
								}
								otherlv_4=':'
								{
									newLeafNode(otherlv_4, grammarAccess.getShapeLayoutAccess().getColonKeyword_2_0_1());
								}
								(
									(
										{
											newCompositeNode(grammarAccess.getShapeLayoutAccess().getXNumberParserRuleCall_2_0_2_0());
										}
										lv_x_5_0=ruleNumber
										{
											if ($current==null) {
												$current = createModelElementForParent(grammarAccess.getShapeLayoutRule());
											}
											set(
												$current,
												"x",
												lv_x_5_0,
												"org.eclipse.elk.graph.text.ElkGraph.Number");
											afterParserOrEnumRuleCall();
										}
									)
								)
								otherlv_6=','
								{
									newLeafNode(otherlv_6, grammarAccess.getShapeLayoutAccess().getCommaKeyword_2_0_3());
								}
								(
									(
										{
											newCompositeNode(grammarAccess.getShapeLayoutAccess().getYNumberParserRuleCall_2_0_4_0());
										}
										lv_y_7_0=ruleNumber
										{
											if ($current==null) {
												$current = createModelElementForParent(grammarAccess.getShapeLayoutRule());
											}
											set(
												$current,
												"y",
												lv_y_7_0,
												"org.eclipse.elk.graph.text.ElkGraph.Number");
											afterParserOrEnumRuleCall();
										}
									)
								)
								))
					{ 
						getUnorderedGroupHelper().returnFromSelection(grammarAccess.getShapeLayoutAccess().getUnorderedGroup_2());
					}
				)
			)|
			(
				{getUnorderedGroupHelper().canSelect(grammarAccess.getShapeLayoutAccess().getUnorderedGroup_2(), 1)}?=>(
					{
						getUnorderedGroupHelper().select(grammarAccess.getShapeLayoutAccess().getUnorderedGroup_2(), 1);
					}
								({true}?=>(otherlv_8='size'
								{
									newLeafNode(otherlv_8, grammarAccess.getShapeLayoutAccess().getSizeKeyword_2_1_0());
								}
								otherlv_9=':'
								{
									newLeafNode(otherlv_9, grammarAccess.getShapeLayoutAccess().getColonKeyword_2_1_1());
								}
								(
									(
										{
											newCompositeNode(grammarAccess.getShapeLayoutAccess().getWidthNumberParserRuleCall_2_1_2_0());
										}
										lv_width_10_0=ruleNumber
										{
											if ($current==null) {
												$current = createModelElementForParent(grammarAccess.getShapeLayoutRule());
											}
											set(
												$current,
												"width",
												lv_width_10_0,
												"org.eclipse.elk.graph.text.ElkGraph.Number");
											afterParserOrEnumRuleCall();
										}
									)
								)
								otherlv_11=','
								{
									newLeafNode(otherlv_11, grammarAccess.getShapeLayoutAccess().getCommaKeyword_2_1_3());
								}
								(
									(
										{
											newCompositeNode(grammarAccess.getShapeLayoutAccess().getHeightNumberParserRuleCall_2_1_4_0());
										}
										lv_height_12_0=ruleNumber
										{
											if ($current==null) {
												$current = createModelElementForParent(grammarAccess.getShapeLayoutRule());
											}
											set(
												$current,
												"height",
												lv_height_12_0,
												"org.eclipse.elk.graph.text.ElkGraph.Number");
											afterParserOrEnumRuleCall();
										}
									)
								)
								))
					{ 
						getUnorderedGroupHelper().returnFromSelection(grammarAccess.getShapeLayoutAccess().getUnorderedGroup_2());
					}
				)
			)
					)*
				)
			)
				{ 
				  getUnorderedGroupHelper().leave(grammarAccess.getShapeLayoutAccess().getUnorderedGroup_2());
				}
		)
		otherlv_13=']'
		{
			newLeafNode(otherlv_13, grammarAccess.getShapeLayoutAccess().getRightSquareBracketKeyword_3());
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
		otherlv_0='edge'
		{
			newLeafNode(otherlv_0, grammarAccess.getElkEdgeAccess().getEdgeKeyword_0());
		}
		(
			(
				(
					lv_identifier_1_0=RULE_ID
					{
						newLeafNode(lv_identifier_1_0, grammarAccess.getElkEdgeAccess().getIdentifierIDTerminalRuleCall_1_0_0());
					}
					{
						if ($current==null) {
							$current = createModelElement(grammarAccess.getElkEdgeRule());
						}
						setWithLastConsumed(
							$current,
							"identifier",
							lv_identifier_1_0,
							"org.eclipse.xtext.common.Terminals.ID");
					}
				)
			)
			otherlv_2=':'
			{
				newLeafNode(otherlv_2, grammarAccess.getElkEdgeAccess().getColonKeyword_1_1());
			}
		)?
		(
			(
				{
					if ($current==null) {
						$current = createModelElement(grammarAccess.getElkEdgeRule());
					}
				}
				{
					newCompositeNode(grammarAccess.getElkEdgeAccess().getSourcesElkConnectableShapeCrossReference_2_0());
				}
				ruleQualifiedId
				{
					afterParserOrEnumRuleCall();
				}
			)
		)
		(
			otherlv_4=','
			{
				newLeafNode(otherlv_4, grammarAccess.getElkEdgeAccess().getCommaKeyword_3_0());
			}
			(
				(
					{
						if ($current==null) {
							$current = createModelElement(grammarAccess.getElkEdgeRule());
						}
					}
					{
						newCompositeNode(grammarAccess.getElkEdgeAccess().getSourcesElkConnectableShapeCrossReference_3_1_0());
					}
					ruleQualifiedId
					{
						afterParserOrEnumRuleCall();
					}
				)
			)
		)*
		otherlv_6='->'
		{
			newLeafNode(otherlv_6, grammarAccess.getElkEdgeAccess().getHyphenMinusGreaterThanSignKeyword_4());
		}
		(
			(
				{
					if ($current==null) {
						$current = createModelElement(grammarAccess.getElkEdgeRule());
					}
				}
				{
					newCompositeNode(grammarAccess.getElkEdgeAccess().getTargetsElkConnectableShapeCrossReference_5_0());
				}
				ruleQualifiedId
				{
					afterParserOrEnumRuleCall();
				}
			)
		)
		(
			otherlv_8=','
			{
				newLeafNode(otherlv_8, grammarAccess.getElkEdgeAccess().getCommaKeyword_6_0());
			}
			(
				(
					{
						if ($current==null) {
							$current = createModelElement(grammarAccess.getElkEdgeRule());
						}
					}
					{
						newCompositeNode(grammarAccess.getElkEdgeAccess().getTargetsElkConnectableShapeCrossReference_6_1_0());
					}
					ruleQualifiedId
					{
						afterParserOrEnumRuleCall();
					}
				)
			)
		)*
		(
			otherlv_10='{'
			{
				newLeafNode(otherlv_10, grammarAccess.getElkEdgeAccess().getLeftCurlyBracketKeyword_7_0());
			}
			(
				{
					if ($current==null) {
						$current = createModelElement(grammarAccess.getElkEdgeRule());
					}
					newCompositeNode(grammarAccess.getElkEdgeAccess().getEdgeLayoutParserRuleCall_7_1());
				}
				this_EdgeLayout_11=ruleEdgeLayout[$current]
				{
					$current = $this_EdgeLayout_11.current;
					afterParserOrEnumRuleCall();
				}
			)?
			(
				(
					{
						newCompositeNode(grammarAccess.getElkEdgeAccess().getPropertiesPropertyParserRuleCall_7_2_0());
					}
					lv_properties_12_0=ruleProperty
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getElkEdgeRule());
						}
						add(
							$current,
							"properties",
							lv_properties_12_0,
							"org.eclipse.elk.graph.text.ElkGraph.Property");
						afterParserOrEnumRuleCall();
					}
				)
			)*
			(
				(
					{
						newCompositeNode(grammarAccess.getElkEdgeAccess().getLabelsElkLabelParserRuleCall_7_3_0());
					}
					lv_labels_13_0=ruleElkLabel
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getElkEdgeRule());
						}
						add(
							$current,
							"labels",
							lv_labels_13_0,
							"org.eclipse.elk.graph.text.ElkGraph.ElkLabel");
						afterParserOrEnumRuleCall();
					}
				)
			)*
			otherlv_14='}'
			{
				newLeafNode(otherlv_14, grammarAccess.getElkEdgeAccess().getRightCurlyBracketKeyword_7_4());
			}
		)?
	)
;


// Rule EdgeLayout
ruleEdgeLayout[EObject in_current]  returns [EObject current=in_current]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		otherlv_0='layout'
		{
			newLeafNode(otherlv_0, grammarAccess.getEdgeLayoutAccess().getLayoutKeyword_0());
		}
		otherlv_1='['
		{
			newLeafNode(otherlv_1, grammarAccess.getEdgeLayoutAccess().getLeftSquareBracketKeyword_1());
		}
		(
			(
				(
					{
						newCompositeNode(grammarAccess.getEdgeLayoutAccess().getSectionsElkSingleEdgeSectionParserRuleCall_2_0_0());
					}
					lv_sections_2_0=ruleElkSingleEdgeSection
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getEdgeLayoutRule());
						}
						add(
							$current,
							"sections",
							lv_sections_2_0,
							"org.eclipse.elk.graph.text.ElkGraph.ElkSingleEdgeSection");
						afterParserOrEnumRuleCall();
					}
				)
			)
			    |
			(
				(
					{
						newCompositeNode(grammarAccess.getEdgeLayoutAccess().getSectionsElkEdgeSectionParserRuleCall_2_1_0());
					}
					lv_sections_3_0=ruleElkEdgeSection
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getEdgeLayoutRule());
						}
						add(
							$current,
							"sections",
							lv_sections_3_0,
							"org.eclipse.elk.graph.text.ElkGraph.ElkEdgeSection");
						afterParserOrEnumRuleCall();
					}
				)
			)+
		)
		otherlv_4=']'
		{
			newLeafNode(otherlv_4, grammarAccess.getEdgeLayoutAccess().getRightSquareBracketKeyword_3());
		}
	)
;

// Entry rule entryRuleElkSingleEdgeSection
entryRuleElkSingleEdgeSection returns [EObject current=null]:
	{ newCompositeNode(grammarAccess.getElkSingleEdgeSectionRule()); }
	iv_ruleElkSingleEdgeSection=ruleElkSingleEdgeSection
	{ $current=$iv_ruleElkSingleEdgeSection.current; }
	EOF;

// Rule ElkSingleEdgeSection
ruleElkSingleEdgeSection returns [EObject current=null]
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
					grammarAccess.getElkSingleEdgeSectionAccess().getElkEdgeSectionAction_0(),
					$current);
			}
		)
		(
			(
				(
					{ 
					  getUnorderedGroupHelper().enter(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0());
					}
					(
						(
				(
					{getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 0)}?=>(
						{
							getUnorderedGroupHelper().select(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 0);
						}
									({true}?=>(otherlv_2='incoming'
									{
										newLeafNode(otherlv_2, grammarAccess.getElkSingleEdgeSectionAccess().getIncomingKeyword_1_0_0_0());
									}
									otherlv_3=':'
									{
										newLeafNode(otherlv_3, grammarAccess.getElkSingleEdgeSectionAccess().getColonKeyword_1_0_0_1());
									}
									(
										(
											{
												if ($current==null) {
													$current = createModelElement(grammarAccess.getElkSingleEdgeSectionRule());
												}
											}
											{
												newCompositeNode(grammarAccess.getElkSingleEdgeSectionAccess().getIncomingShapeElkConnectableShapeCrossReference_1_0_0_2_0());
											}
											ruleQualifiedId
											{
												afterParserOrEnumRuleCall();
											}
										)
									)
									))
						{ 
							getUnorderedGroupHelper().returnFromSelection(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0());
						}
					)
				)|
				(
					{getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 1)}?=>(
						{
							getUnorderedGroupHelper().select(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 1);
						}
									({true}?=>(otherlv_5='outgoing'
									{
										newLeafNode(otherlv_5, grammarAccess.getElkSingleEdgeSectionAccess().getOutgoingKeyword_1_0_1_0());
									}
									otherlv_6=':'
									{
										newLeafNode(otherlv_6, grammarAccess.getElkSingleEdgeSectionAccess().getColonKeyword_1_0_1_1());
									}
									(
										(
											{
												if ($current==null) {
													$current = createModelElement(grammarAccess.getElkSingleEdgeSectionRule());
												}
											}
											{
												newCompositeNode(grammarAccess.getElkSingleEdgeSectionAccess().getOutgoingShapeElkConnectableShapeCrossReference_1_0_1_2_0());
											}
											ruleQualifiedId
											{
												afterParserOrEnumRuleCall();
											}
										)
									)
									))
						{ 
							getUnorderedGroupHelper().returnFromSelection(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0());
						}
					)
				)|
				(
					{getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 2)}?=>(
						{
							getUnorderedGroupHelper().select(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 2);
						}
									({true}?=>(otherlv_8='start'
									{
										newLeafNode(otherlv_8, grammarAccess.getElkSingleEdgeSectionAccess().getStartKeyword_1_0_2_0());
									}
									otherlv_9=':'
									{
										newLeafNode(otherlv_9, grammarAccess.getElkSingleEdgeSectionAccess().getColonKeyword_1_0_2_1());
									}
									(
										(
											{
												newCompositeNode(grammarAccess.getElkSingleEdgeSectionAccess().getStartXNumberParserRuleCall_1_0_2_2_0());
											}
											lv_startX_10_0=ruleNumber
											{
												if ($current==null) {
													$current = createModelElementForParent(grammarAccess.getElkSingleEdgeSectionRule());
												}
												set(
													$current,
													"startX",
													lv_startX_10_0,
													"org.eclipse.elk.graph.text.ElkGraph.Number");
												afterParserOrEnumRuleCall();
											}
										)
									)
									otherlv_11=','
									{
										newLeafNode(otherlv_11, grammarAccess.getElkSingleEdgeSectionAccess().getCommaKeyword_1_0_2_3());
									}
									(
										(
											{
												newCompositeNode(grammarAccess.getElkSingleEdgeSectionAccess().getStartYNumberParserRuleCall_1_0_2_4_0());
											}
											lv_startY_12_0=ruleNumber
											{
												if ($current==null) {
													$current = createModelElementForParent(grammarAccess.getElkSingleEdgeSectionRule());
												}
												set(
													$current,
													"startY",
													lv_startY_12_0,
													"org.eclipse.elk.graph.text.ElkGraph.Number");
												afterParserOrEnumRuleCall();
											}
										)
									)
									))
						{ 
							getUnorderedGroupHelper().returnFromSelection(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0());
						}
					)
				)|
				(
					{getUnorderedGroupHelper().canSelect(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 3)}?=>(
						{
							getUnorderedGroupHelper().select(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0(), 3);
						}
									({true}?=>(otherlv_13='end'
									{
										newLeafNode(otherlv_13, grammarAccess.getElkSingleEdgeSectionAccess().getEndKeyword_1_0_3_0());
									}
									otherlv_14=':'
									{
										newLeafNode(otherlv_14, grammarAccess.getElkSingleEdgeSectionAccess().getColonKeyword_1_0_3_1());
									}
									(
										(
											{
												newCompositeNode(grammarAccess.getElkSingleEdgeSectionAccess().getEndXNumberParserRuleCall_1_0_3_2_0());
											}
											lv_endX_15_0=ruleNumber
											{
												if ($current==null) {
													$current = createModelElementForParent(grammarAccess.getElkSingleEdgeSectionRule());
												}
												set(
													$current,
													"endX",
													lv_endX_15_0,
													"org.eclipse.elk.graph.text.ElkGraph.Number");
												afterParserOrEnumRuleCall();
											}
										)
									)
									otherlv_16=','
									{
										newLeafNode(otherlv_16, grammarAccess.getElkSingleEdgeSectionAccess().getCommaKeyword_1_0_3_3());
									}
									(
										(
											{
												newCompositeNode(grammarAccess.getElkSingleEdgeSectionAccess().getEndYNumberParserRuleCall_1_0_3_4_0());
											}
											lv_endY_17_0=ruleNumber
											{
												if ($current==null) {
													$current = createModelElementForParent(grammarAccess.getElkSingleEdgeSectionRule());
												}
												set(
													$current,
													"endY",
													lv_endY_17_0,
													"org.eclipse.elk.graph.text.ElkGraph.Number");
												afterParserOrEnumRuleCall();
											}
										)
									)
									))
						{ 
							getUnorderedGroupHelper().returnFromSelection(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0());
						}
					)
				)
						)*
					)
				)
					{ 
					  getUnorderedGroupHelper().leave(grammarAccess.getElkSingleEdgeSectionAccess().getUnorderedGroup_1_0());
					}
			)
			(
				otherlv_18='bends'
				{
					newLeafNode(otherlv_18, grammarAccess.getElkSingleEdgeSectionAccess().getBendsKeyword_1_1_0());
				}
				otherlv_19=':'
				{
					newLeafNode(otherlv_19, grammarAccess.getElkSingleEdgeSectionAccess().getColonKeyword_1_1_1());
				}
				(
					(
						{
							newCompositeNode(grammarAccess.getElkSingleEdgeSectionAccess().getBendPointsElkBendPointParserRuleCall_1_1_2_0());
						}
						lv_bendPoints_20_0=ruleElkBendPoint
						{
							if ($current==null) {
								$current = createModelElementForParent(grammarAccess.getElkSingleEdgeSectionRule());
							}
							add(
								$current,
								"bendPoints",
								lv_bendPoints_20_0,
								"org.eclipse.elk.graph.text.ElkGraph.ElkBendPoint");
							afterParserOrEnumRuleCall();
						}
					)
				)
				(
					otherlv_21='|'
					{
						newLeafNode(otherlv_21, grammarAccess.getElkSingleEdgeSectionAccess().getVerticalLineKeyword_1_1_3_0());
					}
					(
						(
							{
								newCompositeNode(grammarAccess.getElkSingleEdgeSectionAccess().getBendPointsElkBendPointParserRuleCall_1_1_3_1_0());
							}
							lv_bendPoints_22_0=ruleElkBendPoint
							{
								if ($current==null) {
									$current = createModelElementForParent(grammarAccess.getElkSingleEdgeSectionRule());
								}
								add(
									$current,
									"bendPoints",
									lv_bendPoints_22_0,
									"org.eclipse.elk.graph.text.ElkGraph.ElkBendPoint");
								afterParserOrEnumRuleCall();
							}
						)
					)
				)*
			)?
			(
				(
					{
						newCompositeNode(grammarAccess.getElkSingleEdgeSectionAccess().getPropertiesPropertyParserRuleCall_1_2_0());
					}
					lv_properties_23_0=ruleProperty
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getElkSingleEdgeSectionRule());
						}
						add(
							$current,
							"properties",
							lv_properties_23_0,
							"org.eclipse.elk.graph.text.ElkGraph.Property");
						afterParserOrEnumRuleCall();
					}
				)
			)*
		)
	)
;

// Entry rule entryRuleElkEdgeSection
entryRuleElkEdgeSection returns [EObject current=null]:
	{ newCompositeNode(grammarAccess.getElkEdgeSectionRule()); }
	iv_ruleElkEdgeSection=ruleElkEdgeSection
	{ $current=$iv_ruleElkEdgeSection.current; }
	EOF;

// Rule ElkEdgeSection
ruleElkEdgeSection returns [EObject current=null]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		otherlv_0='section'
		{
			newLeafNode(otherlv_0, grammarAccess.getElkEdgeSectionAccess().getSectionKeyword_0());
		}
		(
			(
				lv_identifier_1_0=RULE_ID
				{
					newLeafNode(lv_identifier_1_0, grammarAccess.getElkEdgeSectionAccess().getIdentifierIDTerminalRuleCall_1_0());
				}
				{
					if ($current==null) {
						$current = createModelElement(grammarAccess.getElkEdgeSectionRule());
					}
					setWithLastConsumed(
						$current,
						"identifier",
						lv_identifier_1_0,
						"org.eclipse.xtext.common.Terminals.ID");
				}
			)
		)
		(
			otherlv_2='->'
			{
				newLeafNode(otherlv_2, grammarAccess.getElkEdgeSectionAccess().getHyphenMinusGreaterThanSignKeyword_2_0());
			}
			(
				(
					{
						if ($current==null) {
							$current = createModelElement(grammarAccess.getElkEdgeSectionRule());
						}
					}
					otherlv_3=RULE_ID
					{
						newLeafNode(otherlv_3, grammarAccess.getElkEdgeSectionAccess().getOutgoingSectionsElkEdgeSectionCrossReference_2_1_0());
					}
				)
			)
			(
				otherlv_4=','
				{
					newLeafNode(otherlv_4, grammarAccess.getElkEdgeSectionAccess().getCommaKeyword_2_2_0());
				}
				(
					(
						{
							if ($current==null) {
								$current = createModelElement(grammarAccess.getElkEdgeSectionRule());
							}
						}
						otherlv_5=RULE_ID
						{
							newLeafNode(otherlv_5, grammarAccess.getElkEdgeSectionAccess().getOutgoingSectionsElkEdgeSectionCrossReference_2_2_1_0());
						}
					)
				)
			)*
		)?
		otherlv_6='['
		{
			newLeafNode(otherlv_6, grammarAccess.getElkEdgeSectionAccess().getLeftSquareBracketKeyword_3());
		}
		(
			(
				(
					{ 
					  getUnorderedGroupHelper().enter(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0());
					}
					(
						(
				(
					{getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 0)}?=>(
						{
							getUnorderedGroupHelper().select(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 0);
						}
									({true}?=>(otherlv_8='incoming'
									{
										newLeafNode(otherlv_8, grammarAccess.getElkEdgeSectionAccess().getIncomingKeyword_4_0_0_0());
									}
									otherlv_9=':'
									{
										newLeafNode(otherlv_9, grammarAccess.getElkEdgeSectionAccess().getColonKeyword_4_0_0_1());
									}
									(
										(
											{
												if ($current==null) {
													$current = createModelElement(grammarAccess.getElkEdgeSectionRule());
												}
											}
											{
												newCompositeNode(grammarAccess.getElkEdgeSectionAccess().getIncomingShapeElkConnectableShapeCrossReference_4_0_0_2_0());
											}
											ruleQualifiedId
											{
												afterParserOrEnumRuleCall();
											}
										)
									)
									))
						{ 
							getUnorderedGroupHelper().returnFromSelection(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0());
						}
					)
				)|
				(
					{getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 1)}?=>(
						{
							getUnorderedGroupHelper().select(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 1);
						}
									({true}?=>(otherlv_11='outgoing'
									{
										newLeafNode(otherlv_11, grammarAccess.getElkEdgeSectionAccess().getOutgoingKeyword_4_0_1_0());
									}
									otherlv_12=':'
									{
										newLeafNode(otherlv_12, grammarAccess.getElkEdgeSectionAccess().getColonKeyword_4_0_1_1());
									}
									(
										(
											{
												if ($current==null) {
													$current = createModelElement(grammarAccess.getElkEdgeSectionRule());
												}
											}
											{
												newCompositeNode(grammarAccess.getElkEdgeSectionAccess().getOutgoingShapeElkConnectableShapeCrossReference_4_0_1_2_0());
											}
											ruleQualifiedId
											{
												afterParserOrEnumRuleCall();
											}
										)
									)
									))
						{ 
							getUnorderedGroupHelper().returnFromSelection(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0());
						}
					)
				)|
				(
					{getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 2)}?=>(
						{
							getUnorderedGroupHelper().select(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 2);
						}
									({true}?=>(otherlv_14='start'
									{
										newLeafNode(otherlv_14, grammarAccess.getElkEdgeSectionAccess().getStartKeyword_4_0_2_0());
									}
									otherlv_15=':'
									{
										newLeafNode(otherlv_15, grammarAccess.getElkEdgeSectionAccess().getColonKeyword_4_0_2_1());
									}
									(
										(
											{
												newCompositeNode(grammarAccess.getElkEdgeSectionAccess().getStartXNumberParserRuleCall_4_0_2_2_0());
											}
											lv_startX_16_0=ruleNumber
											{
												if ($current==null) {
													$current = createModelElementForParent(grammarAccess.getElkEdgeSectionRule());
												}
												set(
													$current,
													"startX",
													lv_startX_16_0,
													"org.eclipse.elk.graph.text.ElkGraph.Number");
												afterParserOrEnumRuleCall();
											}
										)
									)
									otherlv_17=','
									{
										newLeafNode(otherlv_17, grammarAccess.getElkEdgeSectionAccess().getCommaKeyword_4_0_2_3());
									}
									(
										(
											{
												newCompositeNode(grammarAccess.getElkEdgeSectionAccess().getStartYNumberParserRuleCall_4_0_2_4_0());
											}
											lv_startY_18_0=ruleNumber
											{
												if ($current==null) {
													$current = createModelElementForParent(grammarAccess.getElkEdgeSectionRule());
												}
												set(
													$current,
													"startY",
													lv_startY_18_0,
													"org.eclipse.elk.graph.text.ElkGraph.Number");
												afterParserOrEnumRuleCall();
											}
										)
									)
									))
						{ 
							getUnorderedGroupHelper().returnFromSelection(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0());
						}
					)
				)|
				(
					{getUnorderedGroupHelper().canSelect(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 3)}?=>(
						{
							getUnorderedGroupHelper().select(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0(), 3);
						}
									({true}?=>(otherlv_19='end'
									{
										newLeafNode(otherlv_19, grammarAccess.getElkEdgeSectionAccess().getEndKeyword_4_0_3_0());
									}
									otherlv_20=':'
									{
										newLeafNode(otherlv_20, grammarAccess.getElkEdgeSectionAccess().getColonKeyword_4_0_3_1());
									}
									(
										(
											{
												newCompositeNode(grammarAccess.getElkEdgeSectionAccess().getEndXNumberParserRuleCall_4_0_3_2_0());
											}
											lv_endX_21_0=ruleNumber
											{
												if ($current==null) {
													$current = createModelElementForParent(grammarAccess.getElkEdgeSectionRule());
												}
												set(
													$current,
													"endX",
													lv_endX_21_0,
													"org.eclipse.elk.graph.text.ElkGraph.Number");
												afterParserOrEnumRuleCall();
											}
										)
									)
									otherlv_22=','
									{
										newLeafNode(otherlv_22, grammarAccess.getElkEdgeSectionAccess().getCommaKeyword_4_0_3_3());
									}
									(
										(
											{
												newCompositeNode(grammarAccess.getElkEdgeSectionAccess().getEndYNumberParserRuleCall_4_0_3_4_0());
											}
											lv_endY_23_0=ruleNumber
											{
												if ($current==null) {
													$current = createModelElementForParent(grammarAccess.getElkEdgeSectionRule());
												}
												set(
													$current,
													"endY",
													lv_endY_23_0,
													"org.eclipse.elk.graph.text.ElkGraph.Number");
												afterParserOrEnumRuleCall();
											}
										)
									)
									))
						{ 
							getUnorderedGroupHelper().returnFromSelection(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0());
						}
					)
				)
						)*
					)
				)
					{ 
					  getUnorderedGroupHelper().leave(grammarAccess.getElkEdgeSectionAccess().getUnorderedGroup_4_0());
					}
			)
			(
				otherlv_24='bends'
				{
					newLeafNode(otherlv_24, grammarAccess.getElkEdgeSectionAccess().getBendsKeyword_4_1_0());
				}
				otherlv_25=':'
				{
					newLeafNode(otherlv_25, grammarAccess.getElkEdgeSectionAccess().getColonKeyword_4_1_1());
				}
				(
					(
						{
							newCompositeNode(grammarAccess.getElkEdgeSectionAccess().getBendPointsElkBendPointParserRuleCall_4_1_2_0());
						}
						lv_bendPoints_26_0=ruleElkBendPoint
						{
							if ($current==null) {
								$current = createModelElementForParent(grammarAccess.getElkEdgeSectionRule());
							}
							add(
								$current,
								"bendPoints",
								lv_bendPoints_26_0,
								"org.eclipse.elk.graph.text.ElkGraph.ElkBendPoint");
							afterParserOrEnumRuleCall();
						}
					)
				)
				(
					otherlv_27='|'
					{
						newLeafNode(otherlv_27, grammarAccess.getElkEdgeSectionAccess().getVerticalLineKeyword_4_1_3_0());
					}
					(
						(
							{
								newCompositeNode(grammarAccess.getElkEdgeSectionAccess().getBendPointsElkBendPointParserRuleCall_4_1_3_1_0());
							}
							lv_bendPoints_28_0=ruleElkBendPoint
							{
								if ($current==null) {
									$current = createModelElementForParent(grammarAccess.getElkEdgeSectionRule());
								}
								add(
									$current,
									"bendPoints",
									lv_bendPoints_28_0,
									"org.eclipse.elk.graph.text.ElkGraph.ElkBendPoint");
								afterParserOrEnumRuleCall();
							}
						)
					)
				)*
			)?
			(
				(
					{
						newCompositeNode(grammarAccess.getElkEdgeSectionAccess().getPropertiesPropertyParserRuleCall_4_2_0());
					}
					lv_properties_29_0=ruleProperty
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getElkEdgeSectionRule());
						}
						add(
							$current,
							"properties",
							lv_properties_29_0,
							"org.eclipse.elk.graph.text.ElkGraph.Property");
						afterParserOrEnumRuleCall();
					}
				)
			)*
		)
		otherlv_30=']'
		{
			newLeafNode(otherlv_30, grammarAccess.getElkEdgeSectionAccess().getRightSquareBracketKeyword_5());
		}
	)
;

// Entry rule entryRuleElkBendPoint
entryRuleElkBendPoint returns [EObject current=null]:
	{ newCompositeNode(grammarAccess.getElkBendPointRule()); }
	iv_ruleElkBendPoint=ruleElkBendPoint
	{ $current=$iv_ruleElkBendPoint.current; }
	EOF;

// Rule ElkBendPoint
ruleElkBendPoint returns [EObject current=null]
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
					newCompositeNode(grammarAccess.getElkBendPointAccess().getXNumberParserRuleCall_0_0());
				}
				lv_x_0_0=ruleNumber
				{
					if ($current==null) {
						$current = createModelElementForParent(grammarAccess.getElkBendPointRule());
					}
					set(
						$current,
						"x",
						lv_x_0_0,
						"org.eclipse.elk.graph.text.ElkGraph.Number");
					afterParserOrEnumRuleCall();
				}
			)
		)
		otherlv_1=','
		{
			newLeafNode(otherlv_1, grammarAccess.getElkBendPointAccess().getCommaKeyword_1());
		}
		(
			(
				{
					newCompositeNode(grammarAccess.getElkBendPointAccess().getYNumberParserRuleCall_2_0());
				}
				lv_y_2_0=ruleNumber
				{
					if ($current==null) {
						$current = createModelElementForParent(grammarAccess.getElkBendPointRule());
					}
					set(
						$current,
						"y",
						lv_y_2_0,
						"org.eclipse.elk.graph.text.ElkGraph.Number");
					afterParserOrEnumRuleCall();
				}
			)
		)
	)
;

// Entry rule entryRuleQualifiedId
entryRuleQualifiedId returns [String current=null]:
	{ newCompositeNode(grammarAccess.getQualifiedIdRule()); }
	iv_ruleQualifiedId=ruleQualifiedId
	{ $current=$iv_ruleQualifiedId.current.getText(); }
	EOF;

// Rule QualifiedId
ruleQualifiedId returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		this_ID_0=RULE_ID
		{
			$current.merge(this_ID_0);
		}
		{
			newLeafNode(this_ID_0, grammarAccess.getQualifiedIdAccess().getIDTerminalRuleCall_0());
		}
		(
			kw='.'
			{
				$current.merge(kw);
				newLeafNode(kw, grammarAccess.getQualifiedIdAccess().getFullStopKeyword_1_0());
			}
			this_ID_2=RULE_ID
			{
				$current.merge(this_ID_2);
			}
			{
				newLeafNode(this_ID_2, grammarAccess.getQualifiedIdAccess().getIDTerminalRuleCall_1_1());
			}
		)*
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
						"org.eclipse.elk.graph.text.ElkGraph.PropertyKey");
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
							"org.eclipse.elk.graph.text.ElkGraph.StringValue");
						afterParserOrEnumRuleCall();
					}
				)
			)
			    |
			(
				(
					{
						newCompositeNode(grammarAccess.getPropertyAccess().getValueQualifiedIdValueParserRuleCall_2_1_0());
					}
					lv_value_3_0=ruleQualifiedIdValue
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getPropertyRule());
						}
						set(
							$current,
							"value",
							lv_value_3_0,
							"org.eclipse.elk.graph.text.ElkGraph.QualifiedIdValue");
						afterParserOrEnumRuleCall();
					}
				)
			)
			    |
			(
				(
					{
						newCompositeNode(grammarAccess.getPropertyAccess().getValueNumberValueParserRuleCall_2_2_0());
					}
					lv_value_4_0=ruleNumberValue
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getPropertyRule());
						}
						set(
							$current,
							"value",
							lv_value_4_0,
							"org.eclipse.elk.graph.text.ElkGraph.NumberValue");
						afterParserOrEnumRuleCall();
					}
				)
			)
			    |
			(
				(
					{
						newCompositeNode(grammarAccess.getPropertyAccess().getValueBooleanValueParserRuleCall_2_3_0());
					}
					lv_value_5_0=ruleBooleanValue
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getPropertyRule());
						}
						set(
							$current,
							"value",
							lv_value_5_0,
							"org.eclipse.elk.graph.text.ElkGraph.BooleanValue");
						afterParserOrEnumRuleCall();
					}
				)
			)
			    |
			otherlv_6='null'
			{
				newLeafNode(otherlv_6, grammarAccess.getPropertyAccess().getNullKeyword_2_4());
			}
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
		this_ID_0=RULE_ID
		{
			$current.merge(this_ID_0);
		}
		{
			newLeafNode(this_ID_0, grammarAccess.getPropertyKeyAccess().getIDTerminalRuleCall_0());
		}
		(
			kw='.'
			{
				$current.merge(kw);
				newLeafNode(kw, grammarAccess.getPropertyKeyAccess().getFullStopKeyword_1_0());
			}
			this_ID_2=RULE_ID
			{
				$current.merge(this_ID_2);
			}
			{
				newLeafNode(this_ID_2, grammarAccess.getPropertyKeyAccess().getIDTerminalRuleCall_1_1());
			}
		)*
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

// Entry rule entryRuleQualifiedIdValue
entryRuleQualifiedIdValue returns [String current=null]:
	{ newCompositeNode(grammarAccess.getQualifiedIdValueRule()); }
	iv_ruleQualifiedIdValue=ruleQualifiedIdValue
	{ $current=$iv_ruleQualifiedIdValue.current.getText(); }
	EOF;

// Rule QualifiedIdValue
ruleQualifiedIdValue returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	{
		newCompositeNode(grammarAccess.getQualifiedIdValueAccess().getQualifiedIdParserRuleCall());
	}
	this_QualifiedId_0=ruleQualifiedId
	{
		$current.merge(this_QualifiedId_0);
	}
	{
		afterParserOrEnumRuleCall();
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

RULE_SIGNED_INT : ('+'|'-')? RULE_INT;

RULE_FLOAT : ('+'|'-')? (RULE_INT '.' RULE_INT|RULE_INT ('.' RULE_INT)? ('e'|'E') ('+'|'-')? RULE_INT);

RULE_ID : '^'? ('a'..'z'|'A'..'Z'|'_') ('a'..'z'|'A'..'Z'|'_'|'0'..'9')*;

fragment RULE_INT : ('0'..'9')+;

RULE_STRING : ('"' ('\\' .|~(('\\'|'"')))* '"'|'\'' ('\\' .|~(('\\'|'\'')))* '\'');

RULE_ML_COMMENT : '/*' ( options {greedy=false;} : . )*'*/';

RULE_SL_COMMENT : '//' ~(('\n'|'\r'))* ('\r'? '\n')?;

RULE_WS : (' '|'\t'|'\r'|'\n')+;

RULE_ANY_OTHER : .;
