/*******************************************************************************
 * Copyright (c) 2018 Kiel University and others.
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
grammar InternalGRandom;

options {
	superClass=AbstractInternalAntlrParser;
}

@lexer::header {
package org.eclipse.elk.core.debug.grandom.parser.antlr.internal;

// Hack: Use our own Lexer superclass by means of import. 
// Currently there is no other way to specify the superclass for the lexer.
import org.eclipse.xtext.parser.antlr.Lexer;
}

@parser::header {
package org.eclipse.elk.core.debug.grandom.parser.antlr.internal;

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
import org.eclipse.elk.core.debug.grandom.services.GRandomGrammarAccess;

}

@parser::members {

 	private GRandomGrammarAccess grammarAccess;

    public InternalGRandomParser(TokenStream input, GRandomGrammarAccess grammarAccess) {
        this(input);
        this.grammarAccess = grammarAccess;
        registerRules(grammarAccess.getGrammar());
    }

    @Override
    protected String getFirstRuleName() {
    	return "RandGraph";
   	}

   	@Override
   	protected GRandomGrammarAccess getGrammarAccess() {
   		return grammarAccess;
   	}

}

@rulecatch {
    catch (RecognitionException re) {
        recover(input,re);
        appendSkippedTokens();
    }
}

// Entry rule entryRuleRandGraph
entryRuleRandGraph returns [EObject current=null]:
	{ newCompositeNode(grammarAccess.getRandGraphRule()); }
	iv_ruleRandGraph=ruleRandGraph
	{ $current=$iv_ruleRandGraph.current; }
	EOF;

// Rule RandGraph
ruleRandGraph returns [EObject current=null]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		(
			{
				newCompositeNode(grammarAccess.getRandGraphAccess().getConfigsConfigurationParserRuleCall_0());
			}
			lv_configs_0_0=ruleConfiguration
			{
				if ($current==null) {
					$current = createModelElementForParent(grammarAccess.getRandGraphRule());
				}
				add(
					$current,
					"configs",
					lv_configs_0_0,
					"org.eclipse.elk.core.debug.grandom.GRandom.Configuration");
				afterParserOrEnumRuleCall();
			}
		)
	)*
;

// Entry rule entryRuleConfiguration
entryRuleConfiguration returns [EObject current=null]:
	{ newCompositeNode(grammarAccess.getConfigurationRule()); }
	iv_ruleConfiguration=ruleConfiguration
	{ $current=$iv_ruleConfiguration.current; }
	EOF;

// Rule Configuration
ruleConfiguration returns [EObject current=null]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		otherlv_0='generate'
		{
			newLeafNode(otherlv_0, grammarAccess.getConfigurationAccess().getGenerateKeyword_0());
		}
		(
			(
				lv_samples_1_0=RULE_INT
				{
					newLeafNode(lv_samples_1_0, grammarAccess.getConfigurationAccess().getSamplesINTTerminalRuleCall_1_0());
				}
				{
					if ($current==null) {
						$current = createModelElement(grammarAccess.getConfigurationRule());
					}
					setWithLastConsumed(
						$current,
						"samples",
						lv_samples_1_0,
						"org.eclipse.xtext.common.Terminals.INT");
				}
			)
		)
		(
			(
				{
					newCompositeNode(grammarAccess.getConfigurationAccess().getFormFormEnumRuleCall_2_0());
				}
				lv_form_2_0=ruleForm
				{
					if ($current==null) {
						$current = createModelElementForParent(grammarAccess.getConfigurationRule());
					}
					set(
						$current,
						"form",
						lv_form_2_0,
						"org.eclipse.elk.core.debug.grandom.GRandom.Form");
					afterParserOrEnumRuleCall();
				}
			)
		)
		(
			otherlv_3='{'
			{
				newLeafNode(otherlv_3, grammarAccess.getConfigurationAccess().getLeftCurlyBracketKeyword_3_0());
			}
			(
				(
					{ 
					  getUnorderedGroupHelper().enter(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1());
					}
					(
						(
				(
					{getUnorderedGroupHelper().canSelect(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 0)}?=>(
						{
							getUnorderedGroupHelper().select(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 0);
						}
									({true}?=>((
										{
											newCompositeNode(grammarAccess.getConfigurationAccess().getNodesNodesParserRuleCall_3_1_0_0());
										}
										lv_nodes_5_0=ruleNodes
										{
											if ($current==null) {
												$current = createModelElementForParent(grammarAccess.getConfigurationRule());
											}
											set(
												$current,
												"nodes",
												lv_nodes_5_0,
												"org.eclipse.elk.core.debug.grandom.GRandom.Nodes");
											afterParserOrEnumRuleCall();
										}
									)
									))
						{ 
							getUnorderedGroupHelper().returnFromSelection(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1());
						}
					)
				)|
				(
					{getUnorderedGroupHelper().canSelect(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 1)}?=>(
						{
							getUnorderedGroupHelper().select(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 1);
						}
									({true}?=>((
										{
											newCompositeNode(grammarAccess.getConfigurationAccess().getEdgesEdgesParserRuleCall_3_1_1_0());
										}
										lv_edges_6_0=ruleEdges
										{
											if ($current==null) {
												$current = createModelElementForParent(grammarAccess.getConfigurationRule());
											}
											set(
												$current,
												"edges",
												lv_edges_6_0,
												"org.eclipse.elk.core.debug.grandom.GRandom.Edges");
											afterParserOrEnumRuleCall();
										}
									)
									))
						{ 
							getUnorderedGroupHelper().returnFromSelection(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1());
						}
					)
				)|
				(
					{getUnorderedGroupHelper().canSelect(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 2)}?=>(
						{
							getUnorderedGroupHelper().select(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 2);
						}
									({true}?=>((
										(
											lv_mW_7_0='maxWidth'
											{
												newLeafNode(lv_mW_7_0, grammarAccess.getConfigurationAccess().getMWMaxWidthKeyword_3_1_2_0_0());
											}
											{
												if ($current==null) {
													$current = createModelElement(grammarAccess.getConfigurationRule());
												}
												setWithLastConsumed($current, "mW", true, "maxWidth");
											}
										)
									)
									otherlv_8='='
									{
										newLeafNode(otherlv_8, grammarAccess.getConfigurationAccess().getEqualsSignKeyword_3_1_2_1());
									}
									(
										(
											{
												newCompositeNode(grammarAccess.getConfigurationAccess().getMaxWidthIntegerParserRuleCall_3_1_2_2_0());
											}
											lv_maxWidth_9_0=ruleInteger
											{
												if ($current==null) {
													$current = createModelElementForParent(grammarAccess.getConfigurationRule());
												}
												set(
													$current,
													"maxWidth",
													lv_maxWidth_9_0,
													"org.eclipse.elk.core.debug.grandom.GRandom.Integer");
												afterParserOrEnumRuleCall();
											}
										)
									)
									))
						{ 
							getUnorderedGroupHelper().returnFromSelection(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1());
						}
					)
				)|
				(
					{getUnorderedGroupHelper().canSelect(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 3)}?=>(
						{
							getUnorderedGroupHelper().select(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 3);
						}
									({true}?=>((
										(
											lv_mD_10_0='maxDegree'
											{
												newLeafNode(lv_mD_10_0, grammarAccess.getConfigurationAccess().getMDMaxDegreeKeyword_3_1_3_0_0());
											}
											{
												if ($current==null) {
													$current = createModelElement(grammarAccess.getConfigurationRule());
												}
												setWithLastConsumed($current, "mD", true, "maxDegree");
											}
										)
									)
									otherlv_11='='
									{
										newLeafNode(otherlv_11, grammarAccess.getConfigurationAccess().getEqualsSignKeyword_3_1_3_1());
									}
									(
										(
											{
												newCompositeNode(grammarAccess.getConfigurationAccess().getMaxDegreeIntegerParserRuleCall_3_1_3_2_0());
											}
											lv_maxDegree_12_0=ruleInteger
											{
												if ($current==null) {
													$current = createModelElementForParent(grammarAccess.getConfigurationRule());
												}
												set(
													$current,
													"maxDegree",
													lv_maxDegree_12_0,
													"org.eclipse.elk.core.debug.grandom.GRandom.Integer");
												afterParserOrEnumRuleCall();
											}
										)
									)
									))
						{ 
							getUnorderedGroupHelper().returnFromSelection(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1());
						}
					)
				)|
				(
					{getUnorderedGroupHelper().canSelect(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 4)}?=>(
						{
							getUnorderedGroupHelper().select(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 4);
						}
									({true}?=>((
										(
											lv_pF_13_0='partitionFraction'
											{
												newLeafNode(lv_pF_13_0, grammarAccess.getConfigurationAccess().getPFPartitionFractionKeyword_3_1_4_0_0());
											}
											{
												if ($current==null) {
													$current = createModelElement(grammarAccess.getConfigurationRule());
												}
												setWithLastConsumed($current, "pF", true, "partitionFraction");
											}
										)
									)
									otherlv_14='='
									{
										newLeafNode(otherlv_14, grammarAccess.getConfigurationAccess().getEqualsSignKeyword_3_1_4_1());
									}
									(
										(
											{
												newCompositeNode(grammarAccess.getConfigurationAccess().getFractionDoubleQuantityParserRuleCall_3_1_4_2_0());
											}
											lv_fraction_15_0=ruleDoubleQuantity
											{
												if ($current==null) {
													$current = createModelElementForParent(grammarAccess.getConfigurationRule());
												}
												set(
													$current,
													"fraction",
													lv_fraction_15_0,
													"org.eclipse.elk.core.debug.grandom.GRandom.DoubleQuantity");
												afterParserOrEnumRuleCall();
											}
										)
									)
									))
						{ 
							getUnorderedGroupHelper().returnFromSelection(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1());
						}
					)
				)|
				(
					{getUnorderedGroupHelper().canSelect(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 5)}?=>(
						{
							getUnorderedGroupHelper().select(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 5);
						}
									({true}?=>((
										{
											newCompositeNode(grammarAccess.getConfigurationAccess().getHierarchyHierarchyParserRuleCall_3_1_5_0());
										}
										lv_hierarchy_16_0=ruleHierarchy
										{
											if ($current==null) {
												$current = createModelElementForParent(grammarAccess.getConfigurationRule());
											}
											set(
												$current,
												"hierarchy",
												lv_hierarchy_16_0,
												"org.eclipse.elk.core.debug.grandom.GRandom.Hierarchy");
											afterParserOrEnumRuleCall();
										}
									)
									))
						{ 
							getUnorderedGroupHelper().returnFromSelection(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1());
						}
					)
				)|
				(
					{getUnorderedGroupHelper().canSelect(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 6)}?=>(
						{
							getUnorderedGroupHelper().select(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 6);
						}
									({true}?=>(otherlv_17='seed'
									{
										newLeafNode(otherlv_17, grammarAccess.getConfigurationAccess().getSeedKeyword_3_1_6_0());
									}
									otherlv_18='='
									{
										newLeafNode(otherlv_18, grammarAccess.getConfigurationAccess().getEqualsSignKeyword_3_1_6_1());
									}
									(
										(
											{
												newCompositeNode(grammarAccess.getConfigurationAccess().getSeedIntegerParserRuleCall_3_1_6_2_0());
											}
											lv_seed_19_0=ruleInteger
											{
												if ($current==null) {
													$current = createModelElementForParent(grammarAccess.getConfigurationRule());
												}
												set(
													$current,
													"seed",
													lv_seed_19_0,
													"org.eclipse.elk.core.debug.grandom.GRandom.Integer");
												afterParserOrEnumRuleCall();
											}
										)
									)
									))
						{ 
							getUnorderedGroupHelper().returnFromSelection(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1());
						}
					)
				)|
				(
					{getUnorderedGroupHelper().canSelect(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 7)}?=>(
						{
							getUnorderedGroupHelper().select(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 7);
						}
									({true}?=>(otherlv_20='format'
									{
										newLeafNode(otherlv_20, grammarAccess.getConfigurationAccess().getFormatKeyword_3_1_7_0());
									}
									otherlv_21='='
									{
										newLeafNode(otherlv_21, grammarAccess.getConfigurationAccess().getEqualsSignKeyword_3_1_7_1());
									}
									(
										(
											{
												newCompositeNode(grammarAccess.getConfigurationAccess().getFormatFormatsEnumRuleCall_3_1_7_2_0());
											}
											lv_format_22_0=ruleFormats
											{
												if ($current==null) {
													$current = createModelElementForParent(grammarAccess.getConfigurationRule());
												}
												set(
													$current,
													"format",
													lv_format_22_0,
													"org.eclipse.elk.core.debug.grandom.GRandom.Formats");
												afterParserOrEnumRuleCall();
											}
										)
									)
									))
						{ 
							getUnorderedGroupHelper().returnFromSelection(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1());
						}
					)
				)|
				(
					{getUnorderedGroupHelper().canSelect(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 8)}?=>(
						{
							getUnorderedGroupHelper().select(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1(), 8);
						}
									({true}?=>(otherlv_23='filename'
									{
										newLeafNode(otherlv_23, grammarAccess.getConfigurationAccess().getFilenameKeyword_3_1_8_0());
									}
									otherlv_24='='
									{
										newLeafNode(otherlv_24, grammarAccess.getConfigurationAccess().getEqualsSignKeyword_3_1_8_1());
									}
									(
										(
											lv_filename_25_0=RULE_STRING
											{
												newLeafNode(lv_filename_25_0, grammarAccess.getConfigurationAccess().getFilenameSTRINGTerminalRuleCall_3_1_8_2_0());
											}
											{
												if ($current==null) {
													$current = createModelElement(grammarAccess.getConfigurationRule());
												}
												setWithLastConsumed(
													$current,
													"filename",
													lv_filename_25_0,
													"org.eclipse.xtext.common.Terminals.STRING");
											}
										)
									)
									))
						{ 
							getUnorderedGroupHelper().returnFromSelection(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1());
						}
					)
				)
						)*
					)
				)
					{ 
					  getUnorderedGroupHelper().leave(grammarAccess.getConfigurationAccess().getUnorderedGroup_3_1());
					}
			)
			otherlv_26='}'
			{
				newLeafNode(otherlv_26, grammarAccess.getConfigurationAccess().getRightCurlyBracketKeyword_3_2());
			}
		)?
	)
;

// Entry rule entryRuleHierarchy
entryRuleHierarchy returns [EObject current=null]:
	{ newCompositeNode(grammarAccess.getHierarchyRule()); }
	iv_ruleHierarchy=ruleHierarchy
	{ $current=$iv_ruleHierarchy.current; }
	EOF;

// Rule Hierarchy
ruleHierarchy returns [EObject current=null]
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
					grammarAccess.getHierarchyAccess().getHierarchyAction_0(),
					$current);
			}
		)
		otherlv_1='hierarchy'
		{
			newLeafNode(otherlv_1, grammarAccess.getHierarchyAccess().getHierarchyKeyword_1());
		}
		(
			otherlv_2='{'
			{
				newLeafNode(otherlv_2, grammarAccess.getHierarchyAccess().getLeftCurlyBracketKeyword_2_0());
			}
			(
				(
					{ 
					  getUnorderedGroupHelper().enter(grammarAccess.getHierarchyAccess().getUnorderedGroup_2_1());
					}
					(
						(
				(
					{getUnorderedGroupHelper().canSelect(grammarAccess.getHierarchyAccess().getUnorderedGroup_2_1(), 0)}?=>(
						{
							getUnorderedGroupHelper().select(grammarAccess.getHierarchyAccess().getUnorderedGroup_2_1(), 0);
						}
									({true}?=>(otherlv_4='levels'
									{
										newLeafNode(otherlv_4, grammarAccess.getHierarchyAccess().getLevelsKeyword_2_1_0_0());
									}
									otherlv_5='='
									{
										newLeafNode(otherlv_5, grammarAccess.getHierarchyAccess().getEqualsSignKeyword_2_1_0_1());
									}
									(
										(
											{
												newCompositeNode(grammarAccess.getHierarchyAccess().getLevelsDoubleQuantityParserRuleCall_2_1_0_2_0());
											}
											lv_levels_6_0=ruleDoubleQuantity
											{
												if ($current==null) {
													$current = createModelElementForParent(grammarAccess.getHierarchyRule());
												}
												set(
													$current,
													"levels",
													lv_levels_6_0,
													"org.eclipse.elk.core.debug.grandom.GRandom.DoubleQuantity");
												afterParserOrEnumRuleCall();
											}
										)
									)
									))
						{ 
							getUnorderedGroupHelper().returnFromSelection(grammarAccess.getHierarchyAccess().getUnorderedGroup_2_1());
						}
					)
				)|
				(
					{getUnorderedGroupHelper().canSelect(grammarAccess.getHierarchyAccess().getUnorderedGroup_2_1(), 1)}?=>(
						{
							getUnorderedGroupHelper().select(grammarAccess.getHierarchyAccess().getUnorderedGroup_2_1(), 1);
						}
									({true}?=>(otherlv_7='cross-hierarchy edges'
									{
										newLeafNode(otherlv_7, grammarAccess.getHierarchyAccess().getCrossHierarchyEdgesKeyword_2_1_1_0());
									}
									otherlv_8='='
									{
										newLeafNode(otherlv_8, grammarAccess.getHierarchyAccess().getEqualsSignKeyword_2_1_1_1());
									}
									(
										(
											{
												newCompositeNode(grammarAccess.getHierarchyAccess().getEdgesDoubleQuantityParserRuleCall_2_1_1_2_0());
											}
											lv_edges_9_0=ruleDoubleQuantity
											{
												if ($current==null) {
													$current = createModelElementForParent(grammarAccess.getHierarchyRule());
												}
												set(
													$current,
													"edges",
													lv_edges_9_0,
													"org.eclipse.elk.core.debug.grandom.GRandom.DoubleQuantity");
												afterParserOrEnumRuleCall();
											}
										)
									)
									))
						{ 
							getUnorderedGroupHelper().returnFromSelection(grammarAccess.getHierarchyAccess().getUnorderedGroup_2_1());
						}
					)
				)|
				(
					{getUnorderedGroupHelper().canSelect(grammarAccess.getHierarchyAccess().getUnorderedGroup_2_1(), 2)}?=>(
						{
							getUnorderedGroupHelper().select(grammarAccess.getHierarchyAccess().getUnorderedGroup_2_1(), 2);
						}
									({true}?=>(otherlv_10='compound nodes'
									{
										newLeafNode(otherlv_10, grammarAccess.getHierarchyAccess().getCompoundNodesKeyword_2_1_2_0());
									}
									otherlv_11='='
									{
										newLeafNode(otherlv_11, grammarAccess.getHierarchyAccess().getEqualsSignKeyword_2_1_2_1());
									}
									(
										(
											{
												newCompositeNode(grammarAccess.getHierarchyAccess().getNumHierarchNodesDoubleQuantityParserRuleCall_2_1_2_2_0());
											}
											lv_numHierarchNodes_12_0=ruleDoubleQuantity
											{
												if ($current==null) {
													$current = createModelElementForParent(grammarAccess.getHierarchyRule());
												}
												set(
													$current,
													"numHierarchNodes",
													lv_numHierarchNodes_12_0,
													"org.eclipse.elk.core.debug.grandom.GRandom.DoubleQuantity");
												afterParserOrEnumRuleCall();
											}
										)
									)
									))
						{ 
							getUnorderedGroupHelper().returnFromSelection(grammarAccess.getHierarchyAccess().getUnorderedGroup_2_1());
						}
					)
				)|
				(
					{getUnorderedGroupHelper().canSelect(grammarAccess.getHierarchyAccess().getUnorderedGroup_2_1(), 3)}?=>(
						{
							getUnorderedGroupHelper().select(grammarAccess.getHierarchyAccess().getUnorderedGroup_2_1(), 3);
						}
									({true}?=>(otherlv_13='cross-hierarchy relative edges'
									{
										newLeafNode(otherlv_13, grammarAccess.getHierarchyAccess().getCrossHierarchyRelativeEdgesKeyword_2_1_3_0());
									}
									otherlv_14='='
									{
										newLeafNode(otherlv_14, grammarAccess.getHierarchyAccess().getEqualsSignKeyword_2_1_3_1());
									}
									(
										(
											{
												newCompositeNode(grammarAccess.getHierarchyAccess().getCrossHierarchRelDoubleQuantityParserRuleCall_2_1_3_2_0());
											}
											lv_crossHierarchRel_15_0=ruleDoubleQuantity
											{
												if ($current==null) {
													$current = createModelElementForParent(grammarAccess.getHierarchyRule());
												}
												set(
													$current,
													"crossHierarchRel",
													lv_crossHierarchRel_15_0,
													"org.eclipse.elk.core.debug.grandom.GRandom.DoubleQuantity");
												afterParserOrEnumRuleCall();
											}
										)
									)
									))
						{ 
							getUnorderedGroupHelper().returnFromSelection(grammarAccess.getHierarchyAccess().getUnorderedGroup_2_1());
						}
					)
				)
						)*
					)
				)
					{ 
					  getUnorderedGroupHelper().leave(grammarAccess.getHierarchyAccess().getUnorderedGroup_2_1());
					}
			)
			otherlv_16='}'
			{
				newLeafNode(otherlv_16, grammarAccess.getHierarchyAccess().getRightCurlyBracketKeyword_2_2());
			}
		)?
	)
;

// Entry rule entryRuleEdges
entryRuleEdges returns [EObject current=null]:
	{ newCompositeNode(grammarAccess.getEdgesRule()); }
	iv_ruleEdges=ruleEdges
	{ $current=$iv_ruleEdges.current; }
	EOF;

// Rule Edges
ruleEdges returns [EObject current=null]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		(
			otherlv_0='edges'
			{
				newLeafNode(otherlv_0, grammarAccess.getEdgesAccess().getEdgesKeyword_0_0());
			}
			(
				(
					(
						lv_density_1_0='density'
						{
							newLeafNode(lv_density_1_0, grammarAccess.getEdgesAccess().getDensityDensityKeyword_0_1_0_0());
						}
						{
							if ($current==null) {
								$current = createModelElement(grammarAccess.getEdgesRule());
							}
							setWithLastConsumed($current, "density", true, "density");
						}
					)
				)
				    |
				(
					(
						lv_total_2_0='total'
						{
							newLeafNode(lv_total_2_0, grammarAccess.getEdgesAccess().getTotalTotalKeyword_0_1_1_0());
						}
						{
							if ($current==null) {
								$current = createModelElement(grammarAccess.getEdgesRule());
							}
							setWithLastConsumed($current, "total", true, "total");
						}
					)
				)
				    |
				(
					(
						lv_relative_3_0='relative'
						{
							newLeafNode(lv_relative_3_0, grammarAccess.getEdgesAccess().getRelativeRelativeKeyword_0_1_2_0());
						}
						{
							if ($current==null) {
								$current = createModelElement(grammarAccess.getEdgesRule());
							}
							setWithLastConsumed($current, "relative", true, "relative");
						}
					)
				)
				    |
				(
					(
						lv_outbound_4_0='outgoing'
						{
							newLeafNode(lv_outbound_4_0, grammarAccess.getEdgesAccess().getOutboundOutgoingKeyword_0_1_3_0());
						}
						{
							if ($current==null) {
								$current = createModelElement(grammarAccess.getEdgesRule());
							}
							setWithLastConsumed($current, "outbound", true, "outgoing");
						}
					)
				)
			)
			otherlv_5='='
			{
				newLeafNode(otherlv_5, grammarAccess.getEdgesAccess().getEqualsSignKeyword_0_2());
			}
			(
				(
					{
						newCompositeNode(grammarAccess.getEdgesAccess().getNEdgesDoubleQuantityParserRuleCall_0_3_0());
					}
					lv_nEdges_6_0=ruleDoubleQuantity
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getEdgesRule());
						}
						set(
							$current,
							"nEdges",
							lv_nEdges_6_0,
							"org.eclipse.elk.core.debug.grandom.GRandom.DoubleQuantity");
						afterParserOrEnumRuleCall();
					}
				)
			)
		)
		(
			otherlv_7='{'
			{
				newLeafNode(otherlv_7, grammarAccess.getEdgesAccess().getLeftCurlyBracketKeyword_1_0());
			}
			(
				(
					{ 
					  getUnorderedGroupHelper().enter(grammarAccess.getEdgesAccess().getUnorderedGroup_1_1());
					}
					(
						(
				(
					{getUnorderedGroupHelper().canSelect(grammarAccess.getEdgesAccess().getUnorderedGroup_1_1(), 0)}?=>(
						{
							getUnorderedGroupHelper().select(grammarAccess.getEdgesAccess().getUnorderedGroup_1_1(), 0);
						}
									({true}?=>((
										lv_labels_9_0='labels'
										{
											newLeafNode(lv_labels_9_0, grammarAccess.getEdgesAccess().getLabelsLabelsKeyword_1_1_0_0());
										}
										{
											if ($current==null) {
												$current = createModelElement(grammarAccess.getEdgesRule());
											}
											setWithLastConsumed($current, "labels", true, "labels");
										}
									)
									))
						{ 
							getUnorderedGroupHelper().returnFromSelection(grammarAccess.getEdgesAccess().getUnorderedGroup_1_1());
						}
					)
				)|
				(
					{getUnorderedGroupHelper().canSelect(grammarAccess.getEdgesAccess().getUnorderedGroup_1_1(), 1)}?=>(
						{
							getUnorderedGroupHelper().select(grammarAccess.getEdgesAccess().getUnorderedGroup_1_1(), 1);
						}
									({true}?=>((
										lv_selfLoops_10_0='self loops'
										{
											newLeafNode(lv_selfLoops_10_0, grammarAccess.getEdgesAccess().getSelfLoopsSelfLoopsKeyword_1_1_1_0());
										}
										{
											if ($current==null) {
												$current = createModelElement(grammarAccess.getEdgesRule());
											}
											setWithLastConsumed($current, "selfLoops", true, "self loops");
										}
									)
									))
						{ 
							getUnorderedGroupHelper().returnFromSelection(grammarAccess.getEdgesAccess().getUnorderedGroup_1_1());
						}
					)
				)
						)*
					)
				)
					{ 
					  getUnorderedGroupHelper().leave(grammarAccess.getEdgesAccess().getUnorderedGroup_1_1());
					}
			)
			otherlv_11='}'
			{
				newLeafNode(otherlv_11, grammarAccess.getEdgesAccess().getRightCurlyBracketKeyword_1_2());
			}
		)?
	)
;

// Entry rule entryRuleNodes
entryRuleNodes returns [EObject current=null]:
	{ newCompositeNode(grammarAccess.getNodesRule()); }
	iv_ruleNodes=ruleNodes
	{ $current=$iv_ruleNodes.current; }
	EOF;

// Rule Nodes
ruleNodes returns [EObject current=null]
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
					grammarAccess.getNodesAccess().getNodesAction_0(),
					$current);
			}
		)
		otherlv_1='nodes'
		{
			newLeafNode(otherlv_1, grammarAccess.getNodesAccess().getNodesKeyword_1());
		}
		otherlv_2='='
		{
			newLeafNode(otherlv_2, grammarAccess.getNodesAccess().getEqualsSignKeyword_2());
		}
		(
			(
				{
					newCompositeNode(grammarAccess.getNodesAccess().getNNodesDoubleQuantityParserRuleCall_3_0());
				}
				lv_nNodes_3_0=ruleDoubleQuantity
				{
					if ($current==null) {
						$current = createModelElementForParent(grammarAccess.getNodesRule());
					}
					set(
						$current,
						"nNodes",
						lv_nNodes_3_0,
						"org.eclipse.elk.core.debug.grandom.GRandom.DoubleQuantity");
					afterParserOrEnumRuleCall();
				}
			)
		)
		(
			otherlv_4='{'
			{
				newLeafNode(otherlv_4, grammarAccess.getNodesAccess().getLeftCurlyBracketKeyword_4_0());
			}
			(
				(
					{ 
					  getUnorderedGroupHelper().enter(grammarAccess.getNodesAccess().getUnorderedGroup_4_1());
					}
					(
						(
				(
					{getUnorderedGroupHelper().canSelect(grammarAccess.getNodesAccess().getUnorderedGroup_4_1(), 0)}?=>(
						{
							getUnorderedGroupHelper().select(grammarAccess.getNodesAccess().getUnorderedGroup_4_1(), 0);
						}
									({true}?=>((
										{
											newCompositeNode(grammarAccess.getNodesAccess().getPortsPortsParserRuleCall_4_1_0_0());
										}
										lv_ports_6_0=rulePorts
										{
											if ($current==null) {
												$current = createModelElementForParent(grammarAccess.getNodesRule());
											}
											set(
												$current,
												"ports",
												lv_ports_6_0,
												"org.eclipse.elk.core.debug.grandom.GRandom.Ports");
											afterParserOrEnumRuleCall();
										}
									)
									))
						{ 
							getUnorderedGroupHelper().returnFromSelection(grammarAccess.getNodesAccess().getUnorderedGroup_4_1());
						}
					)
				)|
				(
					{getUnorderedGroupHelper().canSelect(grammarAccess.getNodesAccess().getUnorderedGroup_4_1(), 1)}?=>(
						{
							getUnorderedGroupHelper().select(grammarAccess.getNodesAccess().getUnorderedGroup_4_1(), 1);
						}
									({true}?=>((
										{
											newCompositeNode(grammarAccess.getNodesAccess().getLabelsLabelsParserRuleCall_4_1_1_0());
										}
										lv_labels_7_0=ruleLabels
										{
											if ($current==null) {
												$current = createModelElementForParent(grammarAccess.getNodesRule());
											}
											set(
												$current,
												"labels",
												true,
												"org.eclipse.elk.core.debug.grandom.GRandom.Labels");
											afterParserOrEnumRuleCall();
										}
									)
									))
						{ 
							getUnorderedGroupHelper().returnFromSelection(grammarAccess.getNodesAccess().getUnorderedGroup_4_1());
						}
					)
				)|
				(
					{getUnorderedGroupHelper().canSelect(grammarAccess.getNodesAccess().getUnorderedGroup_4_1(), 2)}?=>(
						{
							getUnorderedGroupHelper().select(grammarAccess.getNodesAccess().getUnorderedGroup_4_1(), 2);
						}
									({true}?=>((
										{
											newCompositeNode(grammarAccess.getNodesAccess().getSizeSizeParserRuleCall_4_1_2_0());
										}
										lv_size_8_0=ruleSize
										{
											if ($current==null) {
												$current = createModelElementForParent(grammarAccess.getNodesRule());
											}
											set(
												$current,
												"size",
												lv_size_8_0,
												"org.eclipse.elk.core.debug.grandom.GRandom.Size");
											afterParserOrEnumRuleCall();
										}
									)
									))
						{ 
							getUnorderedGroupHelper().returnFromSelection(grammarAccess.getNodesAccess().getUnorderedGroup_4_1());
						}
					)
				)|
				(
					{getUnorderedGroupHelper().canSelect(grammarAccess.getNodesAccess().getUnorderedGroup_4_1(), 3)}?=>(
						{
							getUnorderedGroupHelper().select(grammarAccess.getNodesAccess().getUnorderedGroup_4_1(), 3);
						}
									({true}?=>((
										lv_removeIsolated_9_0='remove isolated'
										{
											newLeafNode(lv_removeIsolated_9_0, grammarAccess.getNodesAccess().getRemoveIsolatedRemoveIsolatedKeyword_4_1_3_0());
										}
										{
											if ($current==null) {
												$current = createModelElement(grammarAccess.getNodesRule());
											}
											setWithLastConsumed($current, "removeIsolated", true, "remove isolated");
										}
									)
									))
						{ 
							getUnorderedGroupHelper().returnFromSelection(grammarAccess.getNodesAccess().getUnorderedGroup_4_1());
						}
					)
				)
						)*
					)
				)
					{ 
					  getUnorderedGroupHelper().leave(grammarAccess.getNodesAccess().getUnorderedGroup_4_1());
					}
			)
			otherlv_10='}'
			{
				newLeafNode(otherlv_10, grammarAccess.getNodesAccess().getRightCurlyBracketKeyword_4_2());
			}
		)?
	)
;

// Entry rule entryRuleSize
entryRuleSize returns [EObject current=null]:
	{ newCompositeNode(grammarAccess.getSizeRule()); }
	iv_ruleSize=ruleSize
	{ $current=$iv_ruleSize.current; }
	EOF;

// Rule Size
ruleSize returns [EObject current=null]
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
					grammarAccess.getSizeAccess().getSizeAction_0(),
					$current);
			}
		)
		(
			otherlv_1='size'
			{
				newLeafNode(otherlv_1, grammarAccess.getSizeAccess().getSizeKeyword_1_0());
			}
			(
				otherlv_2='{'
				{
					newLeafNode(otherlv_2, grammarAccess.getSizeAccess().getLeftCurlyBracketKeyword_1_1_0());
				}
				(
					(
						{ 
						  getUnorderedGroupHelper().enter(grammarAccess.getSizeAccess().getUnorderedGroup_1_1_1());
						}
						(
							(
					(
						{getUnorderedGroupHelper().canSelect(grammarAccess.getSizeAccess().getUnorderedGroup_1_1_1(), 0)}?=>(
							{
								getUnorderedGroupHelper().select(grammarAccess.getSizeAccess().getUnorderedGroup_1_1_1(), 0);
							}
										({true}?=>(otherlv_4='height'
										{
											newLeafNode(otherlv_4, grammarAccess.getSizeAccess().getHeightKeyword_1_1_1_0_0());
										}
										otherlv_5='='
										{
											newLeafNode(otherlv_5, grammarAccess.getSizeAccess().getEqualsSignKeyword_1_1_1_0_1());
										}
										(
											(
												{
													newCompositeNode(grammarAccess.getSizeAccess().getHeightDoubleQuantityParserRuleCall_1_1_1_0_2_0());
												}
												lv_height_6_0=ruleDoubleQuantity
												{
													if ($current==null) {
														$current = createModelElementForParent(grammarAccess.getSizeRule());
													}
													set(
														$current,
														"height",
														lv_height_6_0,
														"org.eclipse.elk.core.debug.grandom.GRandom.DoubleQuantity");
													afterParserOrEnumRuleCall();
												}
											)
										)
										))
							{ 
								getUnorderedGroupHelper().returnFromSelection(grammarAccess.getSizeAccess().getUnorderedGroup_1_1_1());
							}
						)
					)|
					(
						{getUnorderedGroupHelper().canSelect(grammarAccess.getSizeAccess().getUnorderedGroup_1_1_1(), 1)}?=>(
							{
								getUnorderedGroupHelper().select(grammarAccess.getSizeAccess().getUnorderedGroup_1_1_1(), 1);
							}
										({true}?=>(otherlv_7='width'
										{
											newLeafNode(otherlv_7, grammarAccess.getSizeAccess().getWidthKeyword_1_1_1_1_0());
										}
										otherlv_8='='
										{
											newLeafNode(otherlv_8, grammarAccess.getSizeAccess().getEqualsSignKeyword_1_1_1_1_1());
										}
										(
											(
												{
													newCompositeNode(grammarAccess.getSizeAccess().getWidthDoubleQuantityParserRuleCall_1_1_1_1_2_0());
												}
												lv_width_9_0=ruleDoubleQuantity
												{
													if ($current==null) {
														$current = createModelElementForParent(grammarAccess.getSizeRule());
													}
													set(
														$current,
														"width",
														lv_width_9_0,
														"org.eclipse.elk.core.debug.grandom.GRandom.DoubleQuantity");
													afterParserOrEnumRuleCall();
												}
											)
										)
										))
							{ 
								getUnorderedGroupHelper().returnFromSelection(grammarAccess.getSizeAccess().getUnorderedGroup_1_1_1());
							}
						)
					)
							)*
						)
					)
						{ 
						  getUnorderedGroupHelper().leave(grammarAccess.getSizeAccess().getUnorderedGroup_1_1_1());
						}
				)
				otherlv_10='}'
				{
					newLeafNode(otherlv_10, grammarAccess.getSizeAccess().getRightCurlyBracketKeyword_1_1_2());
				}
			)?
		)
	)
;

// Entry rule entryRulePorts
entryRulePorts returns [EObject current=null]:
	{ newCompositeNode(grammarAccess.getPortsRule()); }
	iv_rulePorts=rulePorts
	{ $current=$iv_rulePorts.current; }
	EOF;

// Rule Ports
rulePorts returns [EObject current=null]
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
					grammarAccess.getPortsAccess().getPortsAction_0(),
					$current);
			}
		)
		otherlv_1='ports'
		{
			newLeafNode(otherlv_1, grammarAccess.getPortsAccess().getPortsKeyword_1());
		}
		(
			otherlv_2='{'
			{
				newLeafNode(otherlv_2, grammarAccess.getPortsAccess().getLeftCurlyBracketKeyword_2_0());
			}
			(
				(
					{ 
					  getUnorderedGroupHelper().enter(grammarAccess.getPortsAccess().getUnorderedGroup_2_1());
					}
					(
						(
				(
					{getUnorderedGroupHelper().canSelect(grammarAccess.getPortsAccess().getUnorderedGroup_2_1(), 0)}?=>(
						{
							getUnorderedGroupHelper().select(grammarAccess.getPortsAccess().getUnorderedGroup_2_1(), 0);
						}
									({true}?=>((
										{
											newCompositeNode(grammarAccess.getPortsAccess().getLabelsLabelsParserRuleCall_2_1_0_0());
										}
										lv_labels_4_0=ruleLabels
										{
											if ($current==null) {
												$current = createModelElementForParent(grammarAccess.getPortsRule());
											}
											set(
												$current,
												"labels",
												true,
												"org.eclipse.elk.core.debug.grandom.GRandom.Labels");
											afterParserOrEnumRuleCall();
										}
									)
									))
						{ 
							getUnorderedGroupHelper().returnFromSelection(grammarAccess.getPortsAccess().getUnorderedGroup_2_1());
						}
					)
				)|
				(
					{getUnorderedGroupHelper().canSelect(grammarAccess.getPortsAccess().getUnorderedGroup_2_1(), 1)}?=>(
						{
							getUnorderedGroupHelper().select(grammarAccess.getPortsAccess().getUnorderedGroup_2_1(), 1);
						}
									({true}?=>(otherlv_5='re-use'
									{
										newLeafNode(otherlv_5, grammarAccess.getPortsAccess().getReUseKeyword_2_1_1_0());
									}
									otherlv_6='='
									{
										newLeafNode(otherlv_6, grammarAccess.getPortsAccess().getEqualsSignKeyword_2_1_1_1());
									}
									(
										(
											{
												newCompositeNode(grammarAccess.getPortsAccess().getReUseDoubleQuantityParserRuleCall_2_1_1_2_0());
											}
											lv_reUse_7_0=ruleDoubleQuantity
											{
												if ($current==null) {
													$current = createModelElementForParent(grammarAccess.getPortsRule());
												}
												set(
													$current,
													"reUse",
													lv_reUse_7_0,
													"org.eclipse.elk.core.debug.grandom.GRandom.DoubleQuantity");
												afterParserOrEnumRuleCall();
											}
										)
									)
									))
						{ 
							getUnorderedGroupHelper().returnFromSelection(grammarAccess.getPortsAccess().getUnorderedGroup_2_1());
						}
					)
				)|
				(
					{getUnorderedGroupHelper().canSelect(grammarAccess.getPortsAccess().getUnorderedGroup_2_1(), 2)}?=>(
						{
							getUnorderedGroupHelper().select(grammarAccess.getPortsAccess().getUnorderedGroup_2_1(), 2);
						}
									({true}?=>((
										{
											newCompositeNode(grammarAccess.getPortsAccess().getSizeSizeParserRuleCall_2_1_2_0());
										}
										lv_size_8_0=ruleSize
										{
											if ($current==null) {
												$current = createModelElementForParent(grammarAccess.getPortsRule());
											}
											set(
												$current,
												"size",
												lv_size_8_0,
												"org.eclipse.elk.core.debug.grandom.GRandom.Size");
											afterParserOrEnumRuleCall();
										}
									)
									))
						{ 
							getUnorderedGroupHelper().returnFromSelection(grammarAccess.getPortsAccess().getUnorderedGroup_2_1());
						}
					)
				)|
				(
					{getUnorderedGroupHelper().canSelect(grammarAccess.getPortsAccess().getUnorderedGroup_2_1(), 3)}?=>(
						{
							getUnorderedGroupHelper().select(grammarAccess.getPortsAccess().getUnorderedGroup_2_1(), 3);
						}
									({true}?=>(otherlv_9='constraint'
									{
										newLeafNode(otherlv_9, grammarAccess.getPortsAccess().getConstraintKeyword_2_1_3_0());
									}
									otherlv_10='='
									{
										newLeafNode(otherlv_10, grammarAccess.getPortsAccess().getEqualsSignKeyword_2_1_3_1());
									}
									(
										(
											{
												newCompositeNode(grammarAccess.getPortsAccess().getConstraintConstraintTypeEnumRuleCall_2_1_3_2_0());
											}
											lv_constraint_11_0=ruleConstraintType
											{
												if ($current==null) {
													$current = createModelElementForParent(grammarAccess.getPortsRule());
												}
												set(
													$current,
													"constraint",
													lv_constraint_11_0,
													"org.eclipse.elk.core.debug.grandom.GRandom.ConstraintType");
												afterParserOrEnumRuleCall();
											}
										)
									)
									))
						{ 
							getUnorderedGroupHelper().returnFromSelection(grammarAccess.getPortsAccess().getUnorderedGroup_2_1());
						}
					)
				)|
				(
					{getUnorderedGroupHelper().canSelect(grammarAccess.getPortsAccess().getUnorderedGroup_2_1(), 4)}?=>(
						{
							getUnorderedGroupHelper().select(grammarAccess.getPortsAccess().getUnorderedGroup_2_1(), 4);
						}
									({true}?=>((
										{
											newCompositeNode(grammarAccess.getPortsAccess().getFlowFlowParserRuleCall_2_1_4_0());
										}
										lv_flow_12_0=ruleFlow
										{
											if ($current==null) {
												$current = createModelElementForParent(grammarAccess.getPortsRule());
											}
											add(
												$current,
												"flow",
												lv_flow_12_0,
												"org.eclipse.elk.core.debug.grandom.GRandom.Flow");
											afterParserOrEnumRuleCall();
										}
									)
									))+
						{ 
							getUnorderedGroupHelper().returnFromSelection(grammarAccess.getPortsAccess().getUnorderedGroup_2_1());
						}
					)
				)
						)*
					)
				)
					{ 
					  getUnorderedGroupHelper().leave(grammarAccess.getPortsAccess().getUnorderedGroup_2_1());
					}
			)
			otherlv_13='}'
			{
				newLeafNode(otherlv_13, grammarAccess.getPortsAccess().getRightCurlyBracketKeyword_2_2());
			}
		)?
	)
;

// Entry rule entryRuleFlow
entryRuleFlow returns [EObject current=null]:
	{ newCompositeNode(grammarAccess.getFlowRule()); }
	iv_ruleFlow=ruleFlow
	{ $current=$iv_ruleFlow.current; }
	EOF;

// Rule Flow
ruleFlow returns [EObject current=null]
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
					newCompositeNode(grammarAccess.getFlowAccess().getFlowTypeFlowTypeEnumRuleCall_0_0());
				}
				lv_flowType_0_0=ruleFlowType
				{
					if ($current==null) {
						$current = createModelElementForParent(grammarAccess.getFlowRule());
					}
					set(
						$current,
						"flowType",
						lv_flowType_0_0,
						"org.eclipse.elk.core.debug.grandom.GRandom.FlowType");
					afterParserOrEnumRuleCall();
				}
			)
		)
		(
			(
				{
					newCompositeNode(grammarAccess.getFlowAccess().getSideSideEnumRuleCall_1_0());
				}
				lv_side_1_0=ruleSide
				{
					if ($current==null) {
						$current = createModelElementForParent(grammarAccess.getFlowRule());
					}
					set(
						$current,
						"side",
						lv_side_1_0,
						"org.eclipse.elk.core.debug.grandom.GRandom.Side");
					afterParserOrEnumRuleCall();
				}
			)
		)
		otherlv_2='='
		{
			newLeafNode(otherlv_2, grammarAccess.getFlowAccess().getEqualsSignKeyword_2());
		}
		(
			(
				{
					newCompositeNode(grammarAccess.getFlowAccess().getAmountDoubleQuantityParserRuleCall_3_0());
				}
				lv_amount_3_0=ruleDoubleQuantity
				{
					if ($current==null) {
						$current = createModelElementForParent(grammarAccess.getFlowRule());
					}
					set(
						$current,
						"amount",
						lv_amount_3_0,
						"org.eclipse.elk.core.debug.grandom.GRandom.DoubleQuantity");
					afterParserOrEnumRuleCall();
				}
			)
		)
	)
;

// Entry rule entryRuleLabels
entryRuleLabels returns [String current=null]:
	{ newCompositeNode(grammarAccess.getLabelsRule()); }
	iv_ruleLabels=ruleLabels
	{ $current=$iv_ruleLabels.current.getText(); }
	EOF;

// Rule Labels
ruleLabels returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	kw='labels'
	{
		$current.merge(kw);
		newLeafNode(kw, grammarAccess.getLabelsAccess().getLabelsKeyword());
	}
;

// Entry rule entryRuleDoubleQuantity
entryRuleDoubleQuantity returns [EObject current=null]:
	{ newCompositeNode(grammarAccess.getDoubleQuantityRule()); }
	iv_ruleDoubleQuantity=ruleDoubleQuantity
	{ $current=$iv_ruleDoubleQuantity.current; }
	EOF;

// Rule DoubleQuantity
ruleDoubleQuantity returns [EObject current=null]
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
					newCompositeNode(grammarAccess.getDoubleQuantityAccess().getQuantDoubleParserRuleCall_0_0());
				}
				lv_quant_0_0=ruleDouble
				{
					if ($current==null) {
						$current = createModelElementForParent(grammarAccess.getDoubleQuantityRule());
					}
					set(
						$current,
						"quant",
						lv_quant_0_0,
						"org.eclipse.elk.core.debug.grandom.GRandom.Double");
					afterParserOrEnumRuleCall();
				}
			)
		)
		    |
		(
			(
				(
					{
						newCompositeNode(grammarAccess.getDoubleQuantityAccess().getMinDoubleParserRuleCall_1_0_0());
					}
					lv_min_1_0=ruleDouble
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getDoubleQuantityRule());
						}
						set(
							$current,
							"min",
							lv_min_1_0,
							"org.eclipse.elk.core.debug.grandom.GRandom.Double");
						afterParserOrEnumRuleCall();
					}
				)
			)
			(
				(
					lv_minMax_2_0='to'
					{
						newLeafNode(lv_minMax_2_0, grammarAccess.getDoubleQuantityAccess().getMinMaxToKeyword_1_1_0());
					}
					{
						if ($current==null) {
							$current = createModelElement(grammarAccess.getDoubleQuantityRule());
						}
						setWithLastConsumed($current, "minMax", true, "to");
					}
				)
			)
			(
				(
					{
						newCompositeNode(grammarAccess.getDoubleQuantityAccess().getMaxDoubleParserRuleCall_1_2_0());
					}
					lv_max_3_0=ruleDouble
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getDoubleQuantityRule());
						}
						set(
							$current,
							"max",
							lv_max_3_0,
							"org.eclipse.elk.core.debug.grandom.GRandom.Double");
						afterParserOrEnumRuleCall();
					}
				)
			)
		)
		    |
		(
			(
				(
					{
						newCompositeNode(grammarAccess.getDoubleQuantityAccess().getMeanDoubleParserRuleCall_2_0_0());
					}
					lv_mean_4_0=ruleDouble
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getDoubleQuantityRule());
						}
						set(
							$current,
							"mean",
							lv_mean_4_0,
							"org.eclipse.elk.core.debug.grandom.GRandom.Double");
						afterParserOrEnumRuleCall();
					}
				)
			)
			(
				(
					{
						newCompositeNode(grammarAccess.getDoubleQuantityAccess().getGaussianPmParserRuleCall_2_1_0());
					}
					lv_gaussian_5_0=rulePm
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getDoubleQuantityRule());
						}
						set(
							$current,
							"gaussian",
							true,
							"org.eclipse.elk.core.debug.grandom.GRandom.Pm");
						afterParserOrEnumRuleCall();
					}
				)
			)
			(
				(
					{
						newCompositeNode(grammarAccess.getDoubleQuantityAccess().getStddvDoubleParserRuleCall_2_2_0());
					}
					lv_stddv_6_0=ruleDouble
					{
						if ($current==null) {
							$current = createModelElementForParent(grammarAccess.getDoubleQuantityRule());
						}
						set(
							$current,
							"stddv",
							lv_stddv_6_0,
							"org.eclipse.elk.core.debug.grandom.GRandom.Double");
						afterParserOrEnumRuleCall();
					}
				)
			)
		)
	)
;

// Entry rule entryRulePm
entryRulePm returns [String current=null]:
	{ newCompositeNode(grammarAccess.getPmRule()); }
	iv_rulePm=rulePm
	{ $current=$iv_rulePm.current.getText(); }
	EOF;

// Rule Pm
rulePm returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	kw='+/-'
	{
		$current.merge(kw);
		newLeafNode(kw, grammarAccess.getPmAccess().getPlusSignSolidusHyphenMinusKeyword());
	}
;

// Entry rule entryRuleDouble
entryRuleDouble returns [String current=null]:
	{ newCompositeNode(grammarAccess.getDoubleRule()); }
	iv_ruleDouble=ruleDouble
	{ $current=$iv_ruleDouble.current.getText(); }
	EOF;

// Rule Double
ruleDouble returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		this_INT_0=RULE_INT
		{
			$current.merge(this_INT_0);
		}
		{
			newLeafNode(this_INT_0, grammarAccess.getDoubleAccess().getINTTerminalRuleCall_0());
		}
		(
			kw='.'
			{
				$current.merge(kw);
				newLeafNode(kw, grammarAccess.getDoubleAccess().getFullStopKeyword_1_0());
			}
			this_INT_2=RULE_INT
			{
				$current.merge(this_INT_2);
			}
			{
				newLeafNode(this_INT_2, grammarAccess.getDoubleAccess().getINTTerminalRuleCall_1_1());
			}
		)?
	)
;

// Entry rule entryRuleInteger
entryRuleInteger returns [String current=null]:
	{ newCompositeNode(grammarAccess.getIntegerRule()); }
	iv_ruleInteger=ruleInteger
	{ $current=$iv_ruleInteger.current.getText(); }
	EOF;

// Rule Integer
ruleInteger returns [AntlrDatatypeRuleToken current=new AntlrDatatypeRuleToken()]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		this_INT_0=RULE_INT
		{
			$current.merge(this_INT_0);
		}
		{
			newLeafNode(this_INT_0, grammarAccess.getIntegerAccess().getINTTerminalRuleCall_0());
		}
		(
			kw='.'
			{
				$current.merge(kw);
				newLeafNode(kw, grammarAccess.getIntegerAccess().getFullStopKeyword_1_0());
			}
			this_INT_2=RULE_INT
			{
				$current.merge(this_INT_2);
			}
			{
				newLeafNode(this_INT_2, grammarAccess.getIntegerAccess().getINTTerminalRuleCall_1_1());
			}
		)?
	)
;

// Rule Formats
ruleFormats returns [Enumerator current=null]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		(
			enumLiteral_0='elkt'
			{
				$current = grammarAccess.getFormatsAccess().getElktEnumLiteralDeclaration_0().getEnumLiteral().getInstance();
				newLeafNode(enumLiteral_0, grammarAccess.getFormatsAccess().getElktEnumLiteralDeclaration_0());
			}
		)
		    |
		(
			enumLiteral_1='elkg'
			{
				$current = grammarAccess.getFormatsAccess().getElkgEnumLiteralDeclaration_1().getEnumLiteral().getInstance();
				newLeafNode(enumLiteral_1, grammarAccess.getFormatsAccess().getElkgEnumLiteralDeclaration_1());
			}
		)
	)
;

// Rule Form
ruleForm returns [Enumerator current=null]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		(
			enumLiteral_0='trees'
			{
				$current = grammarAccess.getFormAccess().getTreesEnumLiteralDeclaration_0().getEnumLiteral().getInstance();
				newLeafNode(enumLiteral_0, grammarAccess.getFormAccess().getTreesEnumLiteralDeclaration_0());
			}
		)
		    |
		(
			enumLiteral_1='graphs'
			{
				$current = grammarAccess.getFormAccess().getCustomEnumLiteralDeclaration_1().getEnumLiteral().getInstance();
				newLeafNode(enumLiteral_1, grammarAccess.getFormAccess().getCustomEnumLiteralDeclaration_1());
			}
		)
		    |
		(
			enumLiteral_2='bipartite graphs'
			{
				$current = grammarAccess.getFormAccess().getBipartiteEnumLiteralDeclaration_2().getEnumLiteral().getInstance();
				newLeafNode(enumLiteral_2, grammarAccess.getFormAccess().getBipartiteEnumLiteralDeclaration_2());
			}
		)
		    |
		(
			enumLiteral_3='biconnected graphs'
			{
				$current = grammarAccess.getFormAccess().getBiconnectedEnumLiteralDeclaration_3().getEnumLiteral().getInstance();
				newLeafNode(enumLiteral_3, grammarAccess.getFormAccess().getBiconnectedEnumLiteralDeclaration_3());
			}
		)
		    |
		(
			enumLiteral_4='triconnected graphs'
			{
				$current = grammarAccess.getFormAccess().getTriconnectedEnumLiteralDeclaration_4().getEnumLiteral().getInstance();
				newLeafNode(enumLiteral_4, grammarAccess.getFormAccess().getTriconnectedEnumLiteralDeclaration_4());
			}
		)
		    |
		(
			enumLiteral_5='acyclic graphs'
			{
				$current = grammarAccess.getFormAccess().getAcyclicEnumLiteralDeclaration_5().getEnumLiteral().getInstance();
				newLeafNode(enumLiteral_5, grammarAccess.getFormAccess().getAcyclicEnumLiteralDeclaration_5());
			}
		)
	)
;

// Rule Side
ruleSide returns [Enumerator current=null]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		(
			enumLiteral_0='north'
			{
				$current = grammarAccess.getSideAccess().getNorthEnumLiteralDeclaration_0().getEnumLiteral().getInstance();
				newLeafNode(enumLiteral_0, grammarAccess.getSideAccess().getNorthEnumLiteralDeclaration_0());
			}
		)
		    |
		(
			enumLiteral_1='east'
			{
				$current = grammarAccess.getSideAccess().getEastEnumLiteralDeclaration_1().getEnumLiteral().getInstance();
				newLeafNode(enumLiteral_1, grammarAccess.getSideAccess().getEastEnumLiteralDeclaration_1());
			}
		)
		    |
		(
			enumLiteral_2='south'
			{
				$current = grammarAccess.getSideAccess().getSouthEnumLiteralDeclaration_2().getEnumLiteral().getInstance();
				newLeafNode(enumLiteral_2, grammarAccess.getSideAccess().getSouthEnumLiteralDeclaration_2());
			}
		)
		    |
		(
			enumLiteral_3='west'
			{
				$current = grammarAccess.getSideAccess().getWestEnumLiteralDeclaration_3().getEnumLiteral().getInstance();
				newLeafNode(enumLiteral_3, grammarAccess.getSideAccess().getWestEnumLiteralDeclaration_3());
			}
		)
	)
;

// Rule FlowType
ruleFlowType returns [Enumerator current=null]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		(
			enumLiteral_0='incoming'
			{
				$current = grammarAccess.getFlowTypeAccess().getIncomingEnumLiteralDeclaration_0().getEnumLiteral().getInstance();
				newLeafNode(enumLiteral_0, grammarAccess.getFlowTypeAccess().getIncomingEnumLiteralDeclaration_0());
			}
		)
		    |
		(
			enumLiteral_1='outgoing'
			{
				$current = grammarAccess.getFlowTypeAccess().getOutgoingEnumLiteralDeclaration_1().getEnumLiteral().getInstance();
				newLeafNode(enumLiteral_1, grammarAccess.getFlowTypeAccess().getOutgoingEnumLiteralDeclaration_1());
			}
		)
	)
;

// Rule ConstraintType
ruleConstraintType returns [Enumerator current=null]
@init {
	enterRule();
}
@after {
	leaveRule();
}:
	(
		(
			enumLiteral_0='free'
			{
				$current = grammarAccess.getConstraintTypeAccess().getFreeEnumLiteralDeclaration_0().getEnumLiteral().getInstance();
				newLeafNode(enumLiteral_0, grammarAccess.getConstraintTypeAccess().getFreeEnumLiteralDeclaration_0());
			}
		)
		    |
		(
			enumLiteral_1='side'
			{
				$current = grammarAccess.getConstraintTypeAccess().getSideEnumLiteralDeclaration_1().getEnumLiteral().getInstance();
				newLeafNode(enumLiteral_1, grammarAccess.getConstraintTypeAccess().getSideEnumLiteralDeclaration_1());
			}
		)
		    |
		(
			enumLiteral_2='position'
			{
				$current = grammarAccess.getConstraintTypeAccess().getPositionEnumLiteralDeclaration_2().getEnumLiteral().getInstance();
				newLeafNode(enumLiteral_2, grammarAccess.getConstraintTypeAccess().getPositionEnumLiteralDeclaration_2());
			}
		)
		    |
		(
			enumLiteral_3='order'
			{
				$current = grammarAccess.getConstraintTypeAccess().getOrderEnumLiteralDeclaration_3().getEnumLiteral().getInstance();
				newLeafNode(enumLiteral_3, grammarAccess.getConstraintTypeAccess().getOrderEnumLiteralDeclaration_3());
			}
		)
		    |
		(
			enumLiteral_4='ratio'
			{
				$current = grammarAccess.getConstraintTypeAccess().getRatioEnumLiteralDeclaration_4().getEnumLiteral().getInstance();
				newLeafNode(enumLiteral_4, grammarAccess.getConstraintTypeAccess().getRatioEnumLiteralDeclaration_4());
			}
		)
	)
;

RULE_ID : '^'? ('a'..'z'|'A'..'Z'|'_') ('a'..'z'|'A'..'Z'|'_'|'0'..'9')*;

RULE_INT : ('0'..'9')+;

RULE_STRING : ('"' ('\\' .|~(('\\'|'"')))* '"'|'\'' ('\\' .|~(('\\'|'\'')))* '\'');

RULE_ML_COMMENT : '/*' ( options {greedy=false;} : . )*'*/';

RULE_SL_COMMENT : '//' ~(('\n'|'\r'))* ('\r'? '\n')?;

RULE_WS : (' '|'\t'|'\r'|'\n')+;

RULE_ANY_OTHER : .;
