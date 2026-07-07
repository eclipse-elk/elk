/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
grammar InternalGraphvizDot;

options {
	superClass=AbstractInternalAntlrParser;
}

@lexer::header {
package org.eclipse.elk.alg.graphviz.dot.parser.antlr.internal;

// Hack: Use our own Lexer superclass by means of import. 
// Currently there is no other way to specify the superclass for the lexer.
import org.eclipse.xtext.parser.antlr.Lexer;
}

@parser::header {
package org.eclipse.elk.alg.graphviz.dot.parser.antlr.internal;

import org.eclipse.xtext.*;
import org.eclipse.xtext.parser.*;
import org.eclipse.xtext.parser.impl.*;
import org.eclipse.emf.ecore.util.EcoreUtil;
import org.eclipse.emf.ecore.EObject;
import org.eclipse.emf.common.util.Enumerator;
import org.eclipse.xtext.parser.antlr.AbstractInternalAntlrParser;
import org.eclipse.xtext.parser.antlr.XtextTokenStream;
import org.eclipse.xtext.parser.antlr.XtextTokenStream.HiddenTokens;
import org.eclipse.xtext.parser.antlr.AntlrDatatypeRuleToken;
import org.eclipse.elk.alg.graphviz.dot.services.GraphvizDotGrammarAccess;

}

@parser::members {

 	private GraphvizDotGrammarAccess grammarAccess;

    public InternalGraphvizDotParser(TokenStream input, GraphvizDotGrammarAccess grammarAccess) {
        this(input);
        this.grammarAccess = grammarAccess;
        registerRules(grammarAccess.getGrammar());
    }

    @Override
    protected String getFirstRuleName() {
    	return "GraphvizModel";
   	}

   	@Override
   	protected GraphvizDotGrammarAccess getGrammarAccess() {
   		return grammarAccess;
   	}

}

@rulecatch {
    catch (RecognitionException re) {
        recover(input,re);
        appendSkippedTokens();
    }
}

// Entry rule entryRuleGraphvizModel
entryRuleGraphvizModel returns [EObject current=null]:
	{ newCompositeNode(grammarAccess.getGraphvizModelRule()); }
	iv_ruleGraphvizModel=ruleGraphvizModel
	{ $current=$iv_ruleGraphvizModel.current; }
	EOF;

// Rule GraphvizModel
ruleGraphvizModel returns [EObject current=null]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		(
			{
				newCompositeNode(grammarAccess.getGraphvizModelAccess().getGraphsGraphParserRuleCall_0());
			}
			lv_graphs_0_0=ruleGraph
			{
				if ($current==null) {
					$current = createModelElementForParent(grammarAccess.getGraphvizModelRule());
				}
				add(
					$current,
					"graphs",
					lv_graphs_0_0,
					"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.Graph");
				afterParserOrEnumRuleCall();
			}
		)
	)*
;

// Entry rule entryRuleGraph
entryRuleGraph returns [EObject current=null]:
	{ newCompositeNode(grammarAccess.getGraphRule()); }
	iv_ruleGraph=ruleGraph
	{ $current=$iv_ruleGraph.current; }
	EOF;

// Rule Graph
ruleGraph returns [EObject current=null]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		(
			(
				lv_strict_0_0='strict'
				{
					newLeafNode(lv_strict_0_0, grammarAccess.getGraphAccess().getStrictStrictKeyword_0_0());
				}
				{
					if ($current==null) {
						$current = createModelElement(grammarAccess.getGraphRule());
					}
					setWithLastConsumed($current, "strict", true, "strict");
				}
			)
		)?
		(
			(
				{
					newCompositeNode(grammarAccess.getGraphAccess().getTypeGraphTypeEnumRuleCall_1_0());
				}
				lv_type_1_0=ruleGraphType
				{
					if ($current==null) {
						$current = createModelElementForParent(grammarAccess.getGraphRule());
					}
					set(
						$current,
						"type",
						lv_type_1_0,
						"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.GraphType");
					afterParserOrEnumRuleCall();
				}
			)
		)
		(
			(
				{
					newCompositeNode(grammarAccess.getGraphAccess().getNameDotIDParserRuleCall_2_0());
				}
				lv_name_2_0=ruleDotID
				{
					if ($current==null) {
						$current = createModelElementForParent(grammarAccess.getGraphRule());
					}
					set(
						$current,
						"name",
						lv_name_2_0,
						"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.DotID");
					afterParserOrEnumRuleCall();
				}
			)
		)?
		otherlv_3='{'
		{
			newLeafNode(otherlv_3, grammarAccess.getGraphAccess().getLeftCurlyBracketKeyword_3());
		}
		(
			(
				{
					newCompositeNode(grammarAccess.getGraphAccess().getStatementsStatementParserRuleCall_4_0());
				}
				lv_statements_4_0=ruleStatement
				{
					if ($current==null) {
						$current = createModelElementForParent(grammarAccess.getGraphRule());
					}
					add(
						$current,
						"statements",
						lv_statements_4_0,
						"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.Statement");
					afterParserOrEnumRuleCall();
				}
			)
		)*
		otherlv_5='}'
		{
			newLeafNode(otherlv_5, grammarAccess.getGraphAccess().getRightCurlyBracketKeyword_5());
		}
	)
;

// Entry rule entryRuleStatement
entryRuleStatement returns [EObject current=null]:
	{ newCompositeNode(grammarAccess.getStatementRule()); }
	iv_ruleStatement=ruleStatement
	{ $current=$iv_ruleStatement.current; }
	EOF;

// Rule Statement
ruleStatement returns [EObject current=null]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		(
			{
				newCompositeNode(grammarAccess.getStatementAccess().getNodeStatementParserRuleCall_0_0());
			}
			this_NodeStatement_0=ruleNodeStatement
			{
				$current = $this_NodeStatement_0.current;
				afterParserOrEnumRuleCall();
			}
			    |
			{
				newCompositeNode(grammarAccess.getStatementAccess().getEdgeStatementParserRuleCall_0_1());
			}
			this_EdgeStatement_1=ruleEdgeStatement
			{
				$current = $this_EdgeStatement_1.current;
				afterParserOrEnumRuleCall();
			}
			    |
			{
				newCompositeNode(grammarAccess.getStatementAccess().getAttributeStatementParserRuleCall_0_2());
			}
			this_AttributeStatement_2=ruleAttributeStatement
			{
				$current = $this_AttributeStatement_2.current;
				afterParserOrEnumRuleCall();
			}
			    |
			{
				newCompositeNode(grammarAccess.getStatementAccess().getAttributeParserRuleCall_0_3());
			}
			this_Attribute_3=ruleAttribute
			{
				$current = $this_Attribute_3.current;
				afterParserOrEnumRuleCall();
			}
			    |
			{
				newCompositeNode(grammarAccess.getStatementAccess().getSubgraphParserRuleCall_0_4());
			}
			this_Subgraph_4=ruleSubgraph
			{
				$current = $this_Subgraph_4.current;
				afterParserOrEnumRuleCall();
			}
		)
		(
			otherlv_5=';'
			{
				newLeafNode(otherlv_5, grammarAccess.getStatementAccess().getSemicolonKeyword_1());
			}
		)?
	)
;

// Entry rule entryRuleAttribute
entryRuleAttribute returns [EObject current=null]:
	{ newCompositeNode(grammarAccess.getAttributeRule()); }
	iv_ruleAttribute=ruleAttribute
	{ $current=$iv_ruleAttribute.current; }
	EOF;

// Rule Attribute
ruleAttribute returns [EObject current=null]
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
					newCompositeNode(grammarAccess.getAttributeAccess().getNameDotIDParserRuleCall_0_0());
				}
				lv_name_0_0=ruleDotID
				{
					if ($current==null) {
						$current = createModelElementForParent(grammarAccess.getAttributeRule());
					}
					set(
						$current,
						"name",
						lv_name_0_0,
						"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.DotID");
					afterParserOrEnumRuleCall();
				}
			)
		)
		otherlv_1='='
		{
			newLeafNode(otherlv_1, grammarAccess.getAttributeAccess().getEqualsSignKeyword_1());
		}
		(
			(
				{
					newCompositeNode(grammarAccess.getAttributeAccess().getValueDotIDParserRuleCall_2_0());
				}
				lv_value_2_0=ruleDotID
				{
					if ($current==null) {
						$current = createModelElementForParent(grammarAccess.getAttributeRule());
					}
					set(
						$current,
						"value",
						lv_value_2_0,
						"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.DotID");
					afterParserOrEnumRuleCall();
				}
			)
		)
	)
;

// Entry rule entryRuleNodeStatement
entryRuleNodeStatement returns [EObject current=null]:
	{ newCompositeNode(grammarAccess.getNodeStatementRule()); }
	iv_ruleNodeStatement=ruleNodeStatement
	{ $current=$iv_ruleNodeStatement.current; }
	EOF;

// Rule NodeStatement
ruleNodeStatement returns [EObject current=null]
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
					newCompositeNode(grammarAccess.getNodeStatementAccess().getNodeNodeParserRuleCall_0_0());
				}
				lv_node_0_0=ruleNode
				{
					if ($current==null) {
						$current = createModelElementForParent(grammarAccess.getNodeStatementRule());
					}
					set(
						$current,
						"node",
						lv_node_0_0,
						"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.Node");
					afterParserOrEnumRuleCall();
				}
			)
		)
		(
			otherlv_1='['
			{
				newLeafNode(otherlv_1, grammarAccess.getNodeStatementAccess().getLeftSquareBracketKeyword_1_0());
			}
			(
				(
					(
						{
							newCompositeNode(grammarAccess.getNodeStatementAccess().getAttributesListAttributeParserRuleCall_1_1_0_0());
						}
						lv_attributes_2_0=ruleListAttribute
						{
							if ($current==null) {
								$current = createModelElementForParent(grammarAccess.getNodeStatementRule());
							}
							add(
								$current,
								"attributes",
								lv_attributes_2_0,
								"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.ListAttribute");
							afterParserOrEnumRuleCall();
						}
					)
				)
				(
					(
						otherlv_3=','
						{
							newLeafNode(otherlv_3, grammarAccess.getNodeStatementAccess().getCommaKeyword_1_1_1_0());
						}
					)?
					(
						(
							{
								newCompositeNode(grammarAccess.getNodeStatementAccess().getAttributesListAttributeParserRuleCall_1_1_1_1_0());
							}
							lv_attributes_4_0=ruleListAttribute
							{
								if ($current==null) {
									$current = createModelElementForParent(grammarAccess.getNodeStatementRule());
								}
								add(
									$current,
									"attributes",
									lv_attributes_4_0,
									"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.ListAttribute");
								afterParserOrEnumRuleCall();
							}
						)
					)
				)*
			)?
			otherlv_5=']'
			{
				newLeafNode(otherlv_5, grammarAccess.getNodeStatementAccess().getRightSquareBracketKeyword_1_2());
			}
		)?
	)
;

// Entry rule entryRuleNode
entryRuleNode returns [EObject current=null]:
	{ newCompositeNode(grammarAccess.getNodeRule()); }
	iv_ruleNode=ruleNode
	{ $current=$iv_ruleNode.current; }
	EOF;

// Rule Node
ruleNode returns [EObject current=null]
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
					newCompositeNode(grammarAccess.getNodeAccess().getNameDotIDParserRuleCall_0_0());
				}
				lv_name_0_0=ruleDotID
				{
					if ($current==null) {
						$current = createModelElementForParent(grammarAccess.getNodeRule());
					}
					set(
						$current,
						"name",
						lv_name_0_0,
						"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.DotID");
					afterParserOrEnumRuleCall();
				}
			)
		)
		(
			otherlv_1=':'
			{
				newLeafNode(otherlv_1, grammarAccess.getNodeAccess().getColonKeyword_1_0());
			}
			(
				(
					{
						newCompositeNode(grammarAccess.getNodeAccess().getPortPortParserRuleCall_1_1_0());
					}
					lv_port_2_0=rulePort
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getNodeRule());
						}
						set(
							$current,
							"port",
							lv_port_2_0,
							"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.Port");
						afterParserOrEnumRuleCall();
					}
				)
			)
		)?
	)
;

// Entry rule entryRuleEdgeStatement
entryRuleEdgeStatement returns [EObject current=null]:
	{ newCompositeNode(grammarAccess.getEdgeStatementRule()); }
	iv_ruleEdgeStatement=ruleEdgeStatement
	{ $current=$iv_ruleEdgeStatement.current; }
	EOF;

// Rule EdgeStatement
ruleEdgeStatement returns [EObject current=null]
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
					newCompositeNode(grammarAccess.getEdgeStatementAccess().getSourceNodeNodeParserRuleCall_0_0());
				}
				lv_sourceNode_0_0=ruleNode
				{
					if ($current==null) {
						$current = createModelElementForParent(grammarAccess.getEdgeStatementRule());
					}
					set(
						$current,
						"sourceNode",
						lv_sourceNode_0_0,
						"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.Node");
					afterParserOrEnumRuleCall();
				}
			)
		)
		(
			(
				{
					newCompositeNode(grammarAccess.getEdgeStatementAccess().getEdgeTargetsEdgeTargetParserRuleCall_1_0());
				}
				lv_edgeTargets_1_0=ruleEdgeTarget
				{
					if ($current==null) {
						$current = createModelElementForParent(grammarAccess.getEdgeStatementRule());
					}
					add(
						$current,
						"edgeTargets",
						lv_edgeTargets_1_0,
						"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.EdgeTarget");
					afterParserOrEnumRuleCall();
				}
			)
		)+
		(
			otherlv_2='['
			{
				newLeafNode(otherlv_2, grammarAccess.getEdgeStatementAccess().getLeftSquareBracketKeyword_2_0());
			}
			(
				(
					(
						{
							newCompositeNode(grammarAccess.getEdgeStatementAccess().getAttributesListAttributeParserRuleCall_2_1_0_0());
						}
						lv_attributes_3_0=ruleListAttribute
						{
							if ($current==null) {
								$current = createModelElementForParent(grammarAccess.getEdgeStatementRule());
							}
							add(
								$current,
								"attributes",
								lv_attributes_3_0,
								"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.ListAttribute");
							afterParserOrEnumRuleCall();
						}
					)
				)
				(
					(
						otherlv_4=','
						{
							newLeafNode(otherlv_4, grammarAccess.getEdgeStatementAccess().getCommaKeyword_2_1_1_0());
						}
					)?
					(
						(
							{
								newCompositeNode(grammarAccess.getEdgeStatementAccess().getAttributesListAttributeParserRuleCall_2_1_1_1_0());
							}
							lv_attributes_5_0=ruleListAttribute
							{
								if ($current==null) {
									$current = createModelElementForParent(grammarAccess.getEdgeStatementRule());
								}
								add(
									$current,
									"attributes",
									lv_attributes_5_0,
									"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.ListAttribute");
								afterParserOrEnumRuleCall();
							}
						)
					)
				)*
			)?
			otherlv_6=']'
			{
				newLeafNode(otherlv_6, grammarAccess.getEdgeStatementAccess().getRightSquareBracketKeyword_2_2());
			}
		)?
	)
;

// Entry rule entryRuleEdgeTarget
entryRuleEdgeTarget returns [EObject current=null]:
	{ newCompositeNode(grammarAccess.getEdgeTargetRule()); }
	iv_ruleEdgeTarget=ruleEdgeTarget
	{ $current=$iv_ruleEdgeTarget.current; }
	EOF;

// Rule EdgeTarget
ruleEdgeTarget returns [EObject current=null]
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
					newCompositeNode(grammarAccess.getEdgeTargetAccess().getOperatorEdgeOperatorEnumRuleCall_0_0());
				}
				lv_operator_0_0=ruleEdgeOperator
				{
					if ($current==null) {
						$current = createModelElementForParent(grammarAccess.getEdgeTargetRule());
					}
					set(
						$current,
						"operator",
						lv_operator_0_0,
						"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.EdgeOperator");
					afterParserOrEnumRuleCall();
				}
			)
		)
		(
			(
				(
					{
						newCompositeNode(grammarAccess.getEdgeTargetAccess().getTargetSubgraphSubgraphParserRuleCall_1_0_0());
					}
					lv_targetSubgraph_1_0=ruleSubgraph
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getEdgeTargetRule());
						}
						set(
							$current,
							"targetSubgraph",
							lv_targetSubgraph_1_0,
							"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.Subgraph");
						afterParserOrEnumRuleCall();
					}
				)
			)
			    |
			(
				(
					{
						newCompositeNode(grammarAccess.getEdgeTargetAccess().getTargetnodeNodeParserRuleCall_1_1_0());
					}
					lv_targetnode_2_0=ruleNode
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getEdgeTargetRule());
						}
						set(
							$current,
							"targetnode",
							lv_targetnode_2_0,
							"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.Node");
						afterParserOrEnumRuleCall();
					}
				)
			)
		)
	)
;

// Entry rule entryRuleAttributeStatement
entryRuleAttributeStatement returns [EObject current=null]:
	{ newCompositeNode(grammarAccess.getAttributeStatementRule()); }
	iv_ruleAttributeStatement=ruleAttributeStatement
	{ $current=$iv_ruleAttributeStatement.current; }
	EOF;

// Rule AttributeStatement
ruleAttributeStatement returns [EObject current=null]
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
					newCompositeNode(grammarAccess.getAttributeStatementAccess().getTypeAttributeTypeEnumRuleCall_0_0());
				}
				lv_type_0_0=ruleAttributeType
				{
					if ($current==null) {
						$current = createModelElementForParent(grammarAccess.getAttributeStatementRule());
					}
					set(
						$current,
						"type",
						lv_type_0_0,
						"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.AttributeType");
					afterParserOrEnumRuleCall();
				}
			)
		)
		otherlv_1='['
		{
			newLeafNode(otherlv_1, grammarAccess.getAttributeStatementAccess().getLeftSquareBracketKeyword_1());
		}
		(
			(
				(
					{
						newCompositeNode(grammarAccess.getAttributeStatementAccess().getAttributesListAttributeParserRuleCall_2_0_0());
					}
					lv_attributes_2_0=ruleListAttribute
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getAttributeStatementRule());
						}
						add(
							$current,
							"attributes",
							lv_attributes_2_0,
							"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.ListAttribute");
						afterParserOrEnumRuleCall();
					}
				)
			)
			(
				(
					otherlv_3=','
					{
						newLeafNode(otherlv_3, grammarAccess.getAttributeStatementAccess().getCommaKeyword_2_1_0());
					}
				)?
				(
					(
						{
							newCompositeNode(grammarAccess.getAttributeStatementAccess().getAttributesListAttributeParserRuleCall_2_1_1_0());
						}
						lv_attributes_4_0=ruleListAttribute
						{
							if ($current==null) {
								$current = createModelElementForParent(grammarAccess.getAttributeStatementRule());
							}
							add(
								$current,
								"attributes",
								lv_attributes_4_0,
								"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.ListAttribute");
							afterParserOrEnumRuleCall();
						}
					)
				)
			)*
		)?
		otherlv_5=']'
		{
			newLeafNode(otherlv_5, grammarAccess.getAttributeStatementAccess().getRightSquareBracketKeyword_3());
		}
	)
;

// Entry rule entryRuleSubgraph
entryRuleSubgraph returns [EObject current=null]:
	{ newCompositeNode(grammarAccess.getSubgraphRule()); }
	iv_ruleSubgraph=ruleSubgraph
	{ $current=$iv_ruleSubgraph.current; }
	EOF;

// Rule Subgraph
ruleSubgraph returns [EObject current=null]
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
					grammarAccess.getSubgraphAccess().getSubgraphAction_0(),
					$current);
			}
		)
		(
			otherlv_1='subgraph'
			{
				newLeafNode(otherlv_1, grammarAccess.getSubgraphAccess().getSubgraphKeyword_1_0());
			}
			(
				(
					lv_name_2_0=RULE_ID
					{
						newLeafNode(lv_name_2_0, grammarAccess.getSubgraphAccess().getNameIDTerminalRuleCall_1_1_0());
					}
					{
						if ($current==null) {
							$current = createModelElement(grammarAccess.getSubgraphRule());
						}
						setWithLastConsumed(
							$current,
							"name",
							lv_name_2_0,
							"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.ID");
					}
				)
			)?
		)?
		otherlv_3='{'
		{
			newLeafNode(otherlv_3, grammarAccess.getSubgraphAccess().getLeftCurlyBracketKeyword_2());
		}
		(
			(
				{
					newCompositeNode(grammarAccess.getSubgraphAccess().getStatementsStatementParserRuleCall_3_0());
				}
				lv_statements_4_0=ruleStatement
				{
					if ($current==null) {
						$current = createModelElementForParent(grammarAccess.getSubgraphRule());
					}
					add(
						$current,
						"statements",
						lv_statements_4_0,
						"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.Statement");
					afterParserOrEnumRuleCall();
				}
			)
		)*
		otherlv_5='}'
		{
			newLeafNode(otherlv_5, grammarAccess.getSubgraphAccess().getRightCurlyBracketKeyword_4());
		}
	)
;

// Entry rule entryRuleListAttribute
entryRuleListAttribute returns [EObject current=null]:
	{ newCompositeNode(grammarAccess.getListAttributeRule()); }
	iv_ruleListAttribute=ruleListAttribute
	{ $current=$iv_ruleListAttribute.current; }
	EOF;

// Rule ListAttribute
ruleListAttribute returns [EObject current=null]
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
					newCompositeNode(grammarAccess.getListAttributeAccess().getNameDotIDParserRuleCall_0_0());
				}
				lv_name_0_0=ruleDotID
				{
					if ($current==null) {
						$current = createModelElementForParent(grammarAccess.getListAttributeRule());
					}
					set(
						$current,
						"name",
						lv_name_0_0,
						"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.DotID");
					afterParserOrEnumRuleCall();
				}
			)
		)
		(
			otherlv_1='='
			{
				newLeafNode(otherlv_1, grammarAccess.getListAttributeAccess().getEqualsSignKeyword_1_0());
			}
			(
				(
					{
						newCompositeNode(grammarAccess.getListAttributeAccess().getValueDotIDParserRuleCall_1_1_0());
					}
					lv_value_2_0=ruleDotID
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getListAttributeRule());
						}
						set(
							$current,
							"value",
							lv_value_2_0,
							"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.DotID");
						afterParserOrEnumRuleCall();
					}
				)
			)
		)?
	)
;

// Entry rule entryRulePort
entryRulePort returns [EObject current=null]:
	{ newCompositeNode(grammarAccess.getPortRule()); }
	iv_rulePort=rulePort
	{ $current=$iv_rulePort.current; }
	EOF;

// Rule Port
rulePort returns [EObject current=null]
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
					newCompositeNode(grammarAccess.getPortAccess().getNameDotIDParserRuleCall_0_0());
				}
				lv_name_0_0=ruleDotID
				{
					if ($current==null) {
						$current = createModelElementForParent(grammarAccess.getPortRule());
					}
					set(
						$current,
						"name",
						lv_name_0_0,
						"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.DotID");
					afterParserOrEnumRuleCall();
				}
			)
		)
		(
			otherlv_1=':'
			{
				newLeafNode(otherlv_1, grammarAccess.getPortAccess().getColonKeyword_1_0());
			}
			(
				(
					lv_compass_pt_2_0=RULE_ID
					{
						newLeafNode(lv_compass_pt_2_0, grammarAccess.getPortAccess().getCompass_ptIDTerminalRuleCall_1_1_0());
					}
					{
						if ($current==null) {
							$current = createModelElement(grammarAccess.getPortRule());
						}
						setWithLastConsumed(
							$current,
							"compass_pt",
							lv_compass_pt_2_0,
							"org.eclipse.elk.alg.graphviz.dot.GraphvizDot.ID");
					}
				)
			)
		)?
	)
;

// Entry rule entryRuleDotID
entryRuleDotID returns [String current=null]:
	{ newCompositeNode(grammarAccess.getDotIDRule()); }
	iv_ruleDotID=ruleDotID
	{ $current=$iv_ruleDotID.current.getText(); }
	EOF;

// Rule DotID
ruleDotID returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()]
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
			newLeafNode(this_ID_0, grammarAccess.getDotIDAccess().getIDTerminalRuleCall_0());
		}
		    |
		this_INT_1=RULE_INT
		{
			$current.merge(this_INT_1);
		}
		{
			newLeafNode(this_INT_1, grammarAccess.getDotIDAccess().getINTTerminalRuleCall_1());
		}
		    |
		this_FLOAT_2=RULE_FLOAT
		{
			$current.merge(this_FLOAT_2);
		}
		{
			newLeafNode(this_FLOAT_2, grammarAccess.getDotIDAccess().getFLOATTerminalRuleCall_2());
		}
		    |
		this_STRING_3=RULE_STRING
		{
			$current.merge(this_STRING_3);
		}
		{
			newLeafNode(this_STRING_3, grammarAccess.getDotIDAccess().getSTRINGTerminalRuleCall_3());
		}
	)
;

// Rule EdgeOperator
ruleEdgeOperator returns [Enumerator current=null]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		(
			enumLiteral_0='->'
			{
				$current = grammarAccess.getEdgeOperatorAccess().getDirectedEnumLiteralDeclaration_0().getEnumLiteral().getInstance();
				newLeafNode(enumLiteral_0, grammarAccess.getEdgeOperatorAccess().getDirectedEnumLiteralDeclaration_0());
			}
		)
		    |
		(
			enumLiteral_1='--'
			{
				$current = grammarAccess.getEdgeOperatorAccess().getUndirectedEnumLiteralDeclaration_1().getEnumLiteral().getInstance();
				newLeafNode(enumLiteral_1, grammarAccess.getEdgeOperatorAccess().getUndirectedEnumLiteralDeclaration_1());
			}
		)
	)
;

// Rule GraphType
ruleGraphType returns [Enumerator current=null]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		(
			enumLiteral_0='graph'
			{
				$current = grammarAccess.getGraphTypeAccess().getGraphEnumLiteralDeclaration_0().getEnumLiteral().getInstance();
				newLeafNode(enumLiteral_0, grammarAccess.getGraphTypeAccess().getGraphEnumLiteralDeclaration_0());
			}
		)
		    |
		(
			enumLiteral_1='digraph'
			{
				$current = grammarAccess.getGraphTypeAccess().getDigraphEnumLiteralDeclaration_1().getEnumLiteral().getInstance();
				newLeafNode(enumLiteral_1, grammarAccess.getGraphTypeAccess().getDigraphEnumLiteralDeclaration_1());
			}
		)
	)
;

// Rule AttributeType
ruleAttributeType returns [Enumerator current=null]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		(
			enumLiteral_0='graph'
			{
				$current = grammarAccess.getAttributeTypeAccess().getGraphEnumLiteralDeclaration_0().getEnumLiteral().getInstance();
				newLeafNode(enumLiteral_0, grammarAccess.getAttributeTypeAccess().getGraphEnumLiteralDeclaration_0());
			}
		)
		    |
		(
			enumLiteral_1='node'
			{
				$current = grammarAccess.getAttributeTypeAccess().getNodeEnumLiteralDeclaration_1().getEnumLiteral().getInstance();
				newLeafNode(enumLiteral_1, grammarAccess.getAttributeTypeAccess().getNodeEnumLiteralDeclaration_1());
			}
		)
		    |
		(
			enumLiteral_2='edge'
			{
				$current = grammarAccess.getAttributeTypeAccess().getEdgeEnumLiteralDeclaration_2().getEnumLiteral().getInstance();
				newLeafNode(enumLiteral_2, grammarAccess.getAttributeTypeAccess().getEdgeEnumLiteralDeclaration_2());
			}
		)
	)
;

RULE_INT : '-'? ('0'..'9')+;

RULE_ID : ('a'..'z'|'A'..'Z'|'_') ('a'..'z'|'A'..'Z'|'_'|'0'..'9')*;

RULE_FLOAT : '-'? ('0'..'9')* '.' ('0'..'9')+;

RULE_STRING : '"' ('\\' '"'|~('"'))* '"';

RULE_PREC_LINE : '#' ~(('\n'|'\r'))* ('\r'? '\n')?;

RULE_ML_COMMENT : '/*' ( options {greedy=false;} : . )*'*/';

RULE_SL_COMMENT : '//' ~(('\n'|'\r'))* ('\r'? '\n')?;

RULE_WS : (' '|'\t'|'\r'|'\n')+;

RULE_ANY_OTHER : .;
