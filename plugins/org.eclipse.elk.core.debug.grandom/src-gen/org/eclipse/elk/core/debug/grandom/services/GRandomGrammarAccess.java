/*******************************************************************************
 * Copyright (c) 2018 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.core.debug.grandom.services;

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
import org.eclipse.xtext.UnorderedGroup;
import org.eclipse.xtext.common.services.TerminalsGrammarAccess;
import org.eclipse.xtext.service.AbstractElementFinder.AbstractGrammarElementFinder;
import org.eclipse.xtext.service.GrammarProvider;

@Singleton
public class GRandomGrammarAccess extends AbstractGrammarElementFinder {
	
	public class RandGraphElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.core.debug.grandom.GRandom.RandGraph");
		private final Assignment cConfigsAssignment = (Assignment)rule.eContents().get(1);
		private final RuleCall cConfigsConfigurationParserRuleCall_0 = (RuleCall)cConfigsAssignment.eContents().get(0);
		
		//RandGraph:
		//	configs+=Configuration*;
		@Override public ParserRule getRule() { return rule; }
		
		//configs+=Configuration*
		public Assignment getConfigsAssignment() { return cConfigsAssignment; }
		
		//Configuration
		public RuleCall getConfigsConfigurationParserRuleCall_0() { return cConfigsConfigurationParserRuleCall_0; }
	}
	public class ConfigurationElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.core.debug.grandom.GRandom.Configuration");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Keyword cGenerateKeyword_0 = (Keyword)cGroup.eContents().get(0);
		private final Assignment cSamplesAssignment_1 = (Assignment)cGroup.eContents().get(1);
		private final RuleCall cSamplesINTTerminalRuleCall_1_0 = (RuleCall)cSamplesAssignment_1.eContents().get(0);
		private final Assignment cFormAssignment_2 = (Assignment)cGroup.eContents().get(2);
		private final RuleCall cFormFormEnumRuleCall_2_0 = (RuleCall)cFormAssignment_2.eContents().get(0);
		private final Group cGroup_3 = (Group)cGroup.eContents().get(3);
		private final Keyword cLeftCurlyBracketKeyword_3_0 = (Keyword)cGroup_3.eContents().get(0);
		private final UnorderedGroup cUnorderedGroup_3_1 = (UnorderedGroup)cGroup_3.eContents().get(1);
		private final Assignment cNodesAssignment_3_1_0 = (Assignment)cUnorderedGroup_3_1.eContents().get(0);
		private final RuleCall cNodesNodesParserRuleCall_3_1_0_0 = (RuleCall)cNodesAssignment_3_1_0.eContents().get(0);
		private final Assignment cEdgesAssignment_3_1_1 = (Assignment)cUnorderedGroup_3_1.eContents().get(1);
		private final RuleCall cEdgesEdgesParserRuleCall_3_1_1_0 = (RuleCall)cEdgesAssignment_3_1_1.eContents().get(0);
		private final Group cGroup_3_1_2 = (Group)cUnorderedGroup_3_1.eContents().get(2);
		private final Assignment cMWAssignment_3_1_2_0 = (Assignment)cGroup_3_1_2.eContents().get(0);
		private final Keyword cMWMaxWidthKeyword_3_1_2_0_0 = (Keyword)cMWAssignment_3_1_2_0.eContents().get(0);
		private final Keyword cEqualsSignKeyword_3_1_2_1 = (Keyword)cGroup_3_1_2.eContents().get(1);
		private final Assignment cMaxWidthAssignment_3_1_2_2 = (Assignment)cGroup_3_1_2.eContents().get(2);
		private final RuleCall cMaxWidthIntegerParserRuleCall_3_1_2_2_0 = (RuleCall)cMaxWidthAssignment_3_1_2_2.eContents().get(0);
		private final Group cGroup_3_1_3 = (Group)cUnorderedGroup_3_1.eContents().get(3);
		private final Assignment cMDAssignment_3_1_3_0 = (Assignment)cGroup_3_1_3.eContents().get(0);
		private final Keyword cMDMaxDegreeKeyword_3_1_3_0_0 = (Keyword)cMDAssignment_3_1_3_0.eContents().get(0);
		private final Keyword cEqualsSignKeyword_3_1_3_1 = (Keyword)cGroup_3_1_3.eContents().get(1);
		private final Assignment cMaxDegreeAssignment_3_1_3_2 = (Assignment)cGroup_3_1_3.eContents().get(2);
		private final RuleCall cMaxDegreeIntegerParserRuleCall_3_1_3_2_0 = (RuleCall)cMaxDegreeAssignment_3_1_3_2.eContents().get(0);
		private final Group cGroup_3_1_4 = (Group)cUnorderedGroup_3_1.eContents().get(4);
		private final Assignment cPFAssignment_3_1_4_0 = (Assignment)cGroup_3_1_4.eContents().get(0);
		private final Keyword cPFPartitionFractionKeyword_3_1_4_0_0 = (Keyword)cPFAssignment_3_1_4_0.eContents().get(0);
		private final Keyword cEqualsSignKeyword_3_1_4_1 = (Keyword)cGroup_3_1_4.eContents().get(1);
		private final Assignment cFractionAssignment_3_1_4_2 = (Assignment)cGroup_3_1_4.eContents().get(2);
		private final RuleCall cFractionDoubleQuantityParserRuleCall_3_1_4_2_0 = (RuleCall)cFractionAssignment_3_1_4_2.eContents().get(0);
		private final Assignment cHierarchyAssignment_3_1_5 = (Assignment)cUnorderedGroup_3_1.eContents().get(5);
		private final RuleCall cHierarchyHierarchyParserRuleCall_3_1_5_0 = (RuleCall)cHierarchyAssignment_3_1_5.eContents().get(0);
		private final Group cGroup_3_1_6 = (Group)cUnorderedGroup_3_1.eContents().get(6);
		private final Keyword cSeedKeyword_3_1_6_0 = (Keyword)cGroup_3_1_6.eContents().get(0);
		private final Keyword cEqualsSignKeyword_3_1_6_1 = (Keyword)cGroup_3_1_6.eContents().get(1);
		private final Assignment cSeedAssignment_3_1_6_2 = (Assignment)cGroup_3_1_6.eContents().get(2);
		private final RuleCall cSeedIntegerParserRuleCall_3_1_6_2_0 = (RuleCall)cSeedAssignment_3_1_6_2.eContents().get(0);
		private final Group cGroup_3_1_7 = (Group)cUnorderedGroup_3_1.eContents().get(7);
		private final Keyword cFormatKeyword_3_1_7_0 = (Keyword)cGroup_3_1_7.eContents().get(0);
		private final Keyword cEqualsSignKeyword_3_1_7_1 = (Keyword)cGroup_3_1_7.eContents().get(1);
		private final Assignment cFormatAssignment_3_1_7_2 = (Assignment)cGroup_3_1_7.eContents().get(2);
		private final RuleCall cFormatFormatsEnumRuleCall_3_1_7_2_0 = (RuleCall)cFormatAssignment_3_1_7_2.eContents().get(0);
		private final Group cGroup_3_1_8 = (Group)cUnorderedGroup_3_1.eContents().get(8);
		private final Keyword cFilenameKeyword_3_1_8_0 = (Keyword)cGroup_3_1_8.eContents().get(0);
		private final Keyword cEqualsSignKeyword_3_1_8_1 = (Keyword)cGroup_3_1_8.eContents().get(1);
		private final Assignment cFilenameAssignment_3_1_8_2 = (Assignment)cGroup_3_1_8.eContents().get(2);
		private final RuleCall cFilenameSTRINGTerminalRuleCall_3_1_8_2_0 = (RuleCall)cFilenameAssignment_3_1_8_2.eContents().get(0);
		private final Keyword cRightCurlyBracketKeyword_3_2 = (Keyword)cGroup_3.eContents().get(2);
		
		//Configuration:
		//	'generate' samples=INT
		//	form=Form ('{' (nodes=Nodes?
		//	& edges=Edges?
		//	& (mW?='maxWidth' '=' maxWidth=Integer)?
		//	& (mD?='maxDegree' '=' maxDegree=Integer)?
		//	& (pF?='partitionFraction' '=' fraction=DoubleQuantity)?
		//	& hierarchy=Hierarchy?
		//	& ('seed' '=' seed=Integer)?
		//	& ('format' '=' format=Formats)?
		//	& ('filename' '=' filename=STRING)?)
		//	'}')?;
		@Override public ParserRule getRule() { return rule; }
		
		//'generate' samples=INT form=Form ('{' (nodes=Nodes? & edges=Edges? & (mW?='maxWidth' '=' maxWidth=Integer)? &
		//(mD?='maxDegree' '=' maxDegree=Integer)? & (pF?='partitionFraction' '=' fraction=DoubleQuantity)? &
		//hierarchy=Hierarchy? & ('seed' '=' seed=Integer)? & ('format' '=' format=Formats)? & ('filename' '=' filename=STRING)?)
		//'}')?
		public Group getGroup() { return cGroup; }
		
		//'generate'
		public Keyword getGenerateKeyword_0() { return cGenerateKeyword_0; }
		
		//samples=INT
		public Assignment getSamplesAssignment_1() { return cSamplesAssignment_1; }
		
		//INT
		public RuleCall getSamplesINTTerminalRuleCall_1_0() { return cSamplesINTTerminalRuleCall_1_0; }
		
		//form=Form
		public Assignment getFormAssignment_2() { return cFormAssignment_2; }
		
		//Form
		public RuleCall getFormFormEnumRuleCall_2_0() { return cFormFormEnumRuleCall_2_0; }
		
		//('{' (nodes=Nodes? & edges=Edges? & (mW?='maxWidth' '=' maxWidth=Integer)? & (mD?='maxDegree' '=' maxDegree=Integer)? &
		//(pF?='partitionFraction' '=' fraction=DoubleQuantity)? & hierarchy=Hierarchy? & ('seed' '=' seed=Integer)? & ('format'
		//'=' format=Formats)? & ('filename' '=' filename=STRING)?) '}')?
		public Group getGroup_3() { return cGroup_3; }
		
		//'{'
		public Keyword getLeftCurlyBracketKeyword_3_0() { return cLeftCurlyBracketKeyword_3_0; }
		
		//nodes=Nodes? & edges=Edges? & (mW?='maxWidth' '=' maxWidth=Integer)? & (mD?='maxDegree' '=' maxDegree=Integer)? &
		//(pF?='partitionFraction' '=' fraction=DoubleQuantity)? & hierarchy=Hierarchy? & ('seed' '=' seed=Integer)? & ('format'
		//'=' format=Formats)? & ('filename' '=' filename=STRING)?
		public UnorderedGroup getUnorderedGroup_3_1() { return cUnorderedGroup_3_1; }
		
		//nodes=Nodes?
		public Assignment getNodesAssignment_3_1_0() { return cNodesAssignment_3_1_0; }
		
		//Nodes
		public RuleCall getNodesNodesParserRuleCall_3_1_0_0() { return cNodesNodesParserRuleCall_3_1_0_0; }
		
		//edges=Edges?
		public Assignment getEdgesAssignment_3_1_1() { return cEdgesAssignment_3_1_1; }
		
		//Edges
		public RuleCall getEdgesEdgesParserRuleCall_3_1_1_0() { return cEdgesEdgesParserRuleCall_3_1_1_0; }
		
		//(mW?='maxWidth' '=' maxWidth=Integer)?
		public Group getGroup_3_1_2() { return cGroup_3_1_2; }
		
		//mW?='maxWidth'
		public Assignment getMWAssignment_3_1_2_0() { return cMWAssignment_3_1_2_0; }
		
		//'maxWidth'
		public Keyword getMWMaxWidthKeyword_3_1_2_0_0() { return cMWMaxWidthKeyword_3_1_2_0_0; }
		
		//'='
		public Keyword getEqualsSignKeyword_3_1_2_1() { return cEqualsSignKeyword_3_1_2_1; }
		
		//maxWidth=Integer
		public Assignment getMaxWidthAssignment_3_1_2_2() { return cMaxWidthAssignment_3_1_2_2; }
		
		//Integer
		public RuleCall getMaxWidthIntegerParserRuleCall_3_1_2_2_0() { return cMaxWidthIntegerParserRuleCall_3_1_2_2_0; }
		
		//(mD?='maxDegree' '=' maxDegree=Integer)?
		public Group getGroup_3_1_3() { return cGroup_3_1_3; }
		
		//mD?='maxDegree'
		public Assignment getMDAssignment_3_1_3_0() { return cMDAssignment_3_1_3_0; }
		
		//'maxDegree'
		public Keyword getMDMaxDegreeKeyword_3_1_3_0_0() { return cMDMaxDegreeKeyword_3_1_3_0_0; }
		
		//'='
		public Keyword getEqualsSignKeyword_3_1_3_1() { return cEqualsSignKeyword_3_1_3_1; }
		
		//maxDegree=Integer
		public Assignment getMaxDegreeAssignment_3_1_3_2() { return cMaxDegreeAssignment_3_1_3_2; }
		
		//Integer
		public RuleCall getMaxDegreeIntegerParserRuleCall_3_1_3_2_0() { return cMaxDegreeIntegerParserRuleCall_3_1_3_2_0; }
		
		//(pF?='partitionFraction' '=' fraction=DoubleQuantity)?
		public Group getGroup_3_1_4() { return cGroup_3_1_4; }
		
		//pF?='partitionFraction'
		public Assignment getPFAssignment_3_1_4_0() { return cPFAssignment_3_1_4_0; }
		
		//'partitionFraction'
		public Keyword getPFPartitionFractionKeyword_3_1_4_0_0() { return cPFPartitionFractionKeyword_3_1_4_0_0; }
		
		//'='
		public Keyword getEqualsSignKeyword_3_1_4_1() { return cEqualsSignKeyword_3_1_4_1; }
		
		//fraction=DoubleQuantity
		public Assignment getFractionAssignment_3_1_4_2() { return cFractionAssignment_3_1_4_2; }
		
		//DoubleQuantity
		public RuleCall getFractionDoubleQuantityParserRuleCall_3_1_4_2_0() { return cFractionDoubleQuantityParserRuleCall_3_1_4_2_0; }
		
		//hierarchy=Hierarchy?
		public Assignment getHierarchyAssignment_3_1_5() { return cHierarchyAssignment_3_1_5; }
		
		//Hierarchy
		public RuleCall getHierarchyHierarchyParserRuleCall_3_1_5_0() { return cHierarchyHierarchyParserRuleCall_3_1_5_0; }
		
		//('seed' '=' seed=Integer)?
		public Group getGroup_3_1_6() { return cGroup_3_1_6; }
		
		//'seed'
		public Keyword getSeedKeyword_3_1_6_0() { return cSeedKeyword_3_1_6_0; }
		
		//'='
		public Keyword getEqualsSignKeyword_3_1_6_1() { return cEqualsSignKeyword_3_1_6_1; }
		
		//seed=Integer
		public Assignment getSeedAssignment_3_1_6_2() { return cSeedAssignment_3_1_6_2; }
		
		//Integer
		public RuleCall getSeedIntegerParserRuleCall_3_1_6_2_0() { return cSeedIntegerParserRuleCall_3_1_6_2_0; }
		
		//('format' '=' format=Formats)?
		public Group getGroup_3_1_7() { return cGroup_3_1_7; }
		
		//'format'
		public Keyword getFormatKeyword_3_1_7_0() { return cFormatKeyword_3_1_7_0; }
		
		//'='
		public Keyword getEqualsSignKeyword_3_1_7_1() { return cEqualsSignKeyword_3_1_7_1; }
		
		//format=Formats
		public Assignment getFormatAssignment_3_1_7_2() { return cFormatAssignment_3_1_7_2; }
		
		//Formats
		public RuleCall getFormatFormatsEnumRuleCall_3_1_7_2_0() { return cFormatFormatsEnumRuleCall_3_1_7_2_0; }
		
		//('filename' '=' filename=STRING)?
		public Group getGroup_3_1_8() { return cGroup_3_1_8; }
		
		//'filename'
		public Keyword getFilenameKeyword_3_1_8_0() { return cFilenameKeyword_3_1_8_0; }
		
		//'='
		public Keyword getEqualsSignKeyword_3_1_8_1() { return cEqualsSignKeyword_3_1_8_1; }
		
		//filename=STRING
		public Assignment getFilenameAssignment_3_1_8_2() { return cFilenameAssignment_3_1_8_2; }
		
		//STRING
		public RuleCall getFilenameSTRINGTerminalRuleCall_3_1_8_2_0() { return cFilenameSTRINGTerminalRuleCall_3_1_8_2_0; }
		
		//'}'
		public Keyword getRightCurlyBracketKeyword_3_2() { return cRightCurlyBracketKeyword_3_2; }
	}
	public class HierarchyElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.core.debug.grandom.GRandom.Hierarchy");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Action cHierarchyAction_0 = (Action)cGroup.eContents().get(0);
		private final Keyword cHierarchyKeyword_1 = (Keyword)cGroup.eContents().get(1);
		private final Group cGroup_2 = (Group)cGroup.eContents().get(2);
		private final Keyword cLeftCurlyBracketKeyword_2_0 = (Keyword)cGroup_2.eContents().get(0);
		private final UnorderedGroup cUnorderedGroup_2_1 = (UnorderedGroup)cGroup_2.eContents().get(1);
		private final Group cGroup_2_1_0 = (Group)cUnorderedGroup_2_1.eContents().get(0);
		private final Keyword cLevelsKeyword_2_1_0_0 = (Keyword)cGroup_2_1_0.eContents().get(0);
		private final Keyword cEqualsSignKeyword_2_1_0_1 = (Keyword)cGroup_2_1_0.eContents().get(1);
		private final Assignment cLevelsAssignment_2_1_0_2 = (Assignment)cGroup_2_1_0.eContents().get(2);
		private final RuleCall cLevelsDoubleQuantityParserRuleCall_2_1_0_2_0 = (RuleCall)cLevelsAssignment_2_1_0_2.eContents().get(0);
		private final Group cGroup_2_1_1 = (Group)cUnorderedGroup_2_1.eContents().get(1);
		private final Keyword cCrossHierarchyEdgesKeyword_2_1_1_0 = (Keyword)cGroup_2_1_1.eContents().get(0);
		private final Keyword cEqualsSignKeyword_2_1_1_1 = (Keyword)cGroup_2_1_1.eContents().get(1);
		private final Assignment cEdgesAssignment_2_1_1_2 = (Assignment)cGroup_2_1_1.eContents().get(2);
		private final RuleCall cEdgesDoubleQuantityParserRuleCall_2_1_1_2_0 = (RuleCall)cEdgesAssignment_2_1_1_2.eContents().get(0);
		private final Group cGroup_2_1_2 = (Group)cUnorderedGroup_2_1.eContents().get(2);
		private final Keyword cCompoundNodesKeyword_2_1_2_0 = (Keyword)cGroup_2_1_2.eContents().get(0);
		private final Keyword cEqualsSignKeyword_2_1_2_1 = (Keyword)cGroup_2_1_2.eContents().get(1);
		private final Assignment cNumHierarchNodesAssignment_2_1_2_2 = (Assignment)cGroup_2_1_2.eContents().get(2);
		private final RuleCall cNumHierarchNodesDoubleQuantityParserRuleCall_2_1_2_2_0 = (RuleCall)cNumHierarchNodesAssignment_2_1_2_2.eContents().get(0);
		private final Group cGroup_2_1_3 = (Group)cUnorderedGroup_2_1.eContents().get(3);
		private final Keyword cCrossHierarchyRelativeEdgesKeyword_2_1_3_0 = (Keyword)cGroup_2_1_3.eContents().get(0);
		private final Keyword cEqualsSignKeyword_2_1_3_1 = (Keyword)cGroup_2_1_3.eContents().get(1);
		private final Assignment cCrossHierarchRelAssignment_2_1_3_2 = (Assignment)cGroup_2_1_3.eContents().get(2);
		private final RuleCall cCrossHierarchRelDoubleQuantityParserRuleCall_2_1_3_2_0 = (RuleCall)cCrossHierarchRelAssignment_2_1_3_2.eContents().get(0);
		private final Keyword cRightCurlyBracketKeyword_2_2 = (Keyword)cGroup_2.eContents().get(2);
		
		//Hierarchy:
		//	{Hierarchy}
		//	'hierarchy' ('{' (('levels' '=' levels=DoubleQuantity)? // MISSING
		//	& ('cross-hierarchy edges' '=' edges=DoubleQuantity)?
		//	& ('compound nodes' '=' numHierarchNodes=DoubleQuantity)?
		//	& ('cross-hierarchy relative edges' '=' crossHierarchRel=DoubleQuantity)?)
		//	'}')?;
		@Override public ParserRule getRule() { return rule; }
		
		//{Hierarchy} 'hierarchy' ('{' (('levels' '=' levels=DoubleQuantity)? // MISSING
		//& ('cross-hierarchy edges' '=' edges=DoubleQuantity)? & ('compound nodes' '=' numHierarchNodes=DoubleQuantity)? &
		//('cross-hierarchy relative edges' '=' crossHierarchRel=DoubleQuantity)?) '}')?
		public Group getGroup() { return cGroup; }
		
		//{Hierarchy}
		public Action getHierarchyAction_0() { return cHierarchyAction_0; }
		
		//'hierarchy'
		public Keyword getHierarchyKeyword_1() { return cHierarchyKeyword_1; }
		
		//('{' (('levels' '=' levels=DoubleQuantity)? // MISSING
		//& ('cross-hierarchy edges' '=' edges=DoubleQuantity)? & ('compound nodes' '=' numHierarchNodes=DoubleQuantity)? &
		//('cross-hierarchy relative edges' '=' crossHierarchRel=DoubleQuantity)?) '}')?
		public Group getGroup_2() { return cGroup_2; }
		
		//'{'
		public Keyword getLeftCurlyBracketKeyword_2_0() { return cLeftCurlyBracketKeyword_2_0; }
		
		//('levels' '=' levels=DoubleQuantity)? // MISSING
		//& ('cross-hierarchy edges' '=' edges=DoubleQuantity)? & ('compound nodes' '=' numHierarchNodes=DoubleQuantity)? &
		//('cross-hierarchy relative edges' '=' crossHierarchRel=DoubleQuantity)?
		public UnorderedGroup getUnorderedGroup_2_1() { return cUnorderedGroup_2_1; }
		
		//('levels' '=' levels=DoubleQuantity)?
		public Group getGroup_2_1_0() { return cGroup_2_1_0; }
		
		//'levels'
		public Keyword getLevelsKeyword_2_1_0_0() { return cLevelsKeyword_2_1_0_0; }
		
		//'='
		public Keyword getEqualsSignKeyword_2_1_0_1() { return cEqualsSignKeyword_2_1_0_1; }
		
		//levels=DoubleQuantity
		public Assignment getLevelsAssignment_2_1_0_2() { return cLevelsAssignment_2_1_0_2; }
		
		//DoubleQuantity
		public RuleCall getLevelsDoubleQuantityParserRuleCall_2_1_0_2_0() { return cLevelsDoubleQuantityParserRuleCall_2_1_0_2_0; }
		
		//('cross-hierarchy edges' '=' edges=DoubleQuantity)?
		public Group getGroup_2_1_1() { return cGroup_2_1_1; }
		
		//'cross-hierarchy edges'
		public Keyword getCrossHierarchyEdgesKeyword_2_1_1_0() { return cCrossHierarchyEdgesKeyword_2_1_1_0; }
		
		//'='
		public Keyword getEqualsSignKeyword_2_1_1_1() { return cEqualsSignKeyword_2_1_1_1; }
		
		//edges=DoubleQuantity
		public Assignment getEdgesAssignment_2_1_1_2() { return cEdgesAssignment_2_1_1_2; }
		
		//DoubleQuantity
		public RuleCall getEdgesDoubleQuantityParserRuleCall_2_1_1_2_0() { return cEdgesDoubleQuantityParserRuleCall_2_1_1_2_0; }
		
		//('compound nodes' '=' numHierarchNodes=DoubleQuantity)?
		public Group getGroup_2_1_2() { return cGroup_2_1_2; }
		
		//'compound nodes'
		public Keyword getCompoundNodesKeyword_2_1_2_0() { return cCompoundNodesKeyword_2_1_2_0; }
		
		//'='
		public Keyword getEqualsSignKeyword_2_1_2_1() { return cEqualsSignKeyword_2_1_2_1; }
		
		//numHierarchNodes=DoubleQuantity
		public Assignment getNumHierarchNodesAssignment_2_1_2_2() { return cNumHierarchNodesAssignment_2_1_2_2; }
		
		//DoubleQuantity
		public RuleCall getNumHierarchNodesDoubleQuantityParserRuleCall_2_1_2_2_0() { return cNumHierarchNodesDoubleQuantityParserRuleCall_2_1_2_2_0; }
		
		//('cross-hierarchy relative edges' '=' crossHierarchRel=DoubleQuantity)?
		public Group getGroup_2_1_3() { return cGroup_2_1_3; }
		
		//'cross-hierarchy relative edges'
		public Keyword getCrossHierarchyRelativeEdgesKeyword_2_1_3_0() { return cCrossHierarchyRelativeEdgesKeyword_2_1_3_0; }
		
		//'='
		public Keyword getEqualsSignKeyword_2_1_3_1() { return cEqualsSignKeyword_2_1_3_1; }
		
		//crossHierarchRel=DoubleQuantity
		public Assignment getCrossHierarchRelAssignment_2_1_3_2() { return cCrossHierarchRelAssignment_2_1_3_2; }
		
		//DoubleQuantity
		public RuleCall getCrossHierarchRelDoubleQuantityParserRuleCall_2_1_3_2_0() { return cCrossHierarchRelDoubleQuantityParserRuleCall_2_1_3_2_0; }
		
		//'}'
		public Keyword getRightCurlyBracketKeyword_2_2() { return cRightCurlyBracketKeyword_2_2; }
	}
	public class EdgesElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.core.debug.grandom.GRandom.Edges");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Group cGroup_0 = (Group)cGroup.eContents().get(0);
		private final Keyword cEdgesKeyword_0_0 = (Keyword)cGroup_0.eContents().get(0);
		private final Alternatives cAlternatives_0_1 = (Alternatives)cGroup_0.eContents().get(1);
		private final Assignment cDensityAssignment_0_1_0 = (Assignment)cAlternatives_0_1.eContents().get(0);
		private final Keyword cDensityDensityKeyword_0_1_0_0 = (Keyword)cDensityAssignment_0_1_0.eContents().get(0);
		private final Assignment cTotalAssignment_0_1_1 = (Assignment)cAlternatives_0_1.eContents().get(1);
		private final Keyword cTotalTotalKeyword_0_1_1_0 = (Keyword)cTotalAssignment_0_1_1.eContents().get(0);
		private final Assignment cRelativeAssignment_0_1_2 = (Assignment)cAlternatives_0_1.eContents().get(2);
		private final Keyword cRelativeRelativeKeyword_0_1_2_0 = (Keyword)cRelativeAssignment_0_1_2.eContents().get(0);
		private final Assignment cOutboundAssignment_0_1_3 = (Assignment)cAlternatives_0_1.eContents().get(3);
		private final Keyword cOutboundOutgoingKeyword_0_1_3_0 = (Keyword)cOutboundAssignment_0_1_3.eContents().get(0);
		private final Keyword cEqualsSignKeyword_0_2 = (Keyword)cGroup_0.eContents().get(2);
		private final Assignment cNEdgesAssignment_0_3 = (Assignment)cGroup_0.eContents().get(3);
		private final RuleCall cNEdgesDoubleQuantityParserRuleCall_0_3_0 = (RuleCall)cNEdgesAssignment_0_3.eContents().get(0);
		private final Group cGroup_1 = (Group)cGroup.eContents().get(1);
		private final Keyword cLeftCurlyBracketKeyword_1_0 = (Keyword)cGroup_1.eContents().get(0);
		private final UnorderedGroup cUnorderedGroup_1_1 = (UnorderedGroup)cGroup_1.eContents().get(1);
		private final Assignment cLabelsAssignment_1_1_0 = (Assignment)cUnorderedGroup_1_1.eContents().get(0);
		private final Keyword cLabelsLabelsKeyword_1_1_0_0 = (Keyword)cLabelsAssignment_1_1_0.eContents().get(0);
		private final Assignment cSelfLoopsAssignment_1_1_1 = (Assignment)cUnorderedGroup_1_1.eContents().get(1);
		private final Keyword cSelfLoopsSelfLoopsKeyword_1_1_1_0 = (Keyword)cSelfLoopsAssignment_1_1_1.eContents().get(0);
		private final Keyword cRightCurlyBracketKeyword_1_2 = (Keyword)cGroup_1.eContents().get(2);
		
		//Edges:
		//	('edges' (density?='density' | total?='total' | relative?='relative' | outbound?='outgoing')
		//	'=' nEdges=DoubleQuantity) ('{' (labels?='labels'? & selfLoops?='self loops'?)
		//	'}')?;
		@Override public ParserRule getRule() { return rule; }
		
		//('edges' (density?='density' | total?='total' | relative?='relative' | outbound?='outgoing') '=' nEdges=DoubleQuantity)
		//('{' (labels?='labels'? & selfLoops?='self loops'?) '}')?
		public Group getGroup() { return cGroup; }
		
		//'edges' (density?='density' | total?='total' | relative?='relative' | outbound?='outgoing') '=' nEdges=DoubleQuantity
		public Group getGroup_0() { return cGroup_0; }
		
		//'edges'
		public Keyword getEdgesKeyword_0_0() { return cEdgesKeyword_0_0; }
		
		//density?='density' | total?='total' | relative?='relative' | outbound?='outgoing'
		public Alternatives getAlternatives_0_1() { return cAlternatives_0_1; }
		
		//density?='density'
		public Assignment getDensityAssignment_0_1_0() { return cDensityAssignment_0_1_0; }
		
		//'density'
		public Keyword getDensityDensityKeyword_0_1_0_0() { return cDensityDensityKeyword_0_1_0_0; }
		
		//total?='total'
		public Assignment getTotalAssignment_0_1_1() { return cTotalAssignment_0_1_1; }
		
		//'total'
		public Keyword getTotalTotalKeyword_0_1_1_0() { return cTotalTotalKeyword_0_1_1_0; }
		
		//relative?='relative'
		public Assignment getRelativeAssignment_0_1_2() { return cRelativeAssignment_0_1_2; }
		
		//'relative'
		public Keyword getRelativeRelativeKeyword_0_1_2_0() { return cRelativeRelativeKeyword_0_1_2_0; }
		
		//outbound?='outgoing'
		public Assignment getOutboundAssignment_0_1_3() { return cOutboundAssignment_0_1_3; }
		
		//'outgoing'
		public Keyword getOutboundOutgoingKeyword_0_1_3_0() { return cOutboundOutgoingKeyword_0_1_3_0; }
		
		//'='
		public Keyword getEqualsSignKeyword_0_2() { return cEqualsSignKeyword_0_2; }
		
		//nEdges=DoubleQuantity
		public Assignment getNEdgesAssignment_0_3() { return cNEdgesAssignment_0_3; }
		
		//DoubleQuantity
		public RuleCall getNEdgesDoubleQuantityParserRuleCall_0_3_0() { return cNEdgesDoubleQuantityParserRuleCall_0_3_0; }
		
		//('{' (labels?='labels'? & selfLoops?='self loops'?) '}')?
		public Group getGroup_1() { return cGroup_1; }
		
		//'{'
		public Keyword getLeftCurlyBracketKeyword_1_0() { return cLeftCurlyBracketKeyword_1_0; }
		
		//labels?='labels'? & selfLoops?='self loops'?
		public UnorderedGroup getUnorderedGroup_1_1() { return cUnorderedGroup_1_1; }
		
		//labels?='labels'?
		public Assignment getLabelsAssignment_1_1_0() { return cLabelsAssignment_1_1_0; }
		
		//'labels'
		public Keyword getLabelsLabelsKeyword_1_1_0_0() { return cLabelsLabelsKeyword_1_1_0_0; }
		
		//selfLoops?='self loops'?
		public Assignment getSelfLoopsAssignment_1_1_1() { return cSelfLoopsAssignment_1_1_1; }
		
		//'self loops'
		public Keyword getSelfLoopsSelfLoopsKeyword_1_1_1_0() { return cSelfLoopsSelfLoopsKeyword_1_1_1_0; }
		
		//'}'
		public Keyword getRightCurlyBracketKeyword_1_2() { return cRightCurlyBracketKeyword_1_2; }
	}
	public class NodesElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.core.debug.grandom.GRandom.Nodes");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Action cNodesAction_0 = (Action)cGroup.eContents().get(0);
		private final Keyword cNodesKeyword_1 = (Keyword)cGroup.eContents().get(1);
		private final Keyword cEqualsSignKeyword_2 = (Keyword)cGroup.eContents().get(2);
		private final Assignment cNNodesAssignment_3 = (Assignment)cGroup.eContents().get(3);
		private final RuleCall cNNodesDoubleQuantityParserRuleCall_3_0 = (RuleCall)cNNodesAssignment_3.eContents().get(0);
		private final Group cGroup_4 = (Group)cGroup.eContents().get(4);
		private final Keyword cLeftCurlyBracketKeyword_4_0 = (Keyword)cGroup_4.eContents().get(0);
		private final UnorderedGroup cUnorderedGroup_4_1 = (UnorderedGroup)cGroup_4.eContents().get(1);
		private final Assignment cPortsAssignment_4_1_0 = (Assignment)cUnorderedGroup_4_1.eContents().get(0);
		private final RuleCall cPortsPortsParserRuleCall_4_1_0_0 = (RuleCall)cPortsAssignment_4_1_0.eContents().get(0);
		private final Assignment cLabelsAssignment_4_1_1 = (Assignment)cUnorderedGroup_4_1.eContents().get(1);
		private final RuleCall cLabelsLabelsParserRuleCall_4_1_1_0 = (RuleCall)cLabelsAssignment_4_1_1.eContents().get(0);
		private final Assignment cSizeAssignment_4_1_2 = (Assignment)cUnorderedGroup_4_1.eContents().get(2);
		private final RuleCall cSizeSizeParserRuleCall_4_1_2_0 = (RuleCall)cSizeAssignment_4_1_2.eContents().get(0);
		private final Assignment cRemoveIsolatedAssignment_4_1_3 = (Assignment)cUnorderedGroup_4_1.eContents().get(3);
		private final Keyword cRemoveIsolatedRemoveIsolatedKeyword_4_1_3_0 = (Keyword)cRemoveIsolatedAssignment_4_1_3.eContents().get(0);
		private final Keyword cRightCurlyBracketKeyword_4_2 = (Keyword)cGroup_4.eContents().get(2);
		
		//Nodes:
		//	{Nodes}
		//	'nodes' '=' nNodes=DoubleQuantity ('{' (ports=Ports? & labels?=Labels? & size=Size? &
		//	removeIsolated?='remove isolated'?)
		//	'}')?;
		@Override public ParserRule getRule() { return rule; }
		
		//{Nodes} 'nodes' '=' nNodes=DoubleQuantity ('{' (ports=Ports? & labels?=Labels? & size=Size? &
		//removeIsolated?='remove isolated'?) '}')?
		public Group getGroup() { return cGroup; }
		
		//{Nodes}
		public Action getNodesAction_0() { return cNodesAction_0; }
		
		//'nodes'
		public Keyword getNodesKeyword_1() { return cNodesKeyword_1; }
		
		//'='
		public Keyword getEqualsSignKeyword_2() { return cEqualsSignKeyword_2; }
		
		//nNodes=DoubleQuantity
		public Assignment getNNodesAssignment_3() { return cNNodesAssignment_3; }
		
		//DoubleQuantity
		public RuleCall getNNodesDoubleQuantityParserRuleCall_3_0() { return cNNodesDoubleQuantityParserRuleCall_3_0; }
		
		//('{' (ports=Ports? & labels?=Labels? & size=Size? & removeIsolated?='remove isolated'?) '}')?
		public Group getGroup_4() { return cGroup_4; }
		
		//'{'
		public Keyword getLeftCurlyBracketKeyword_4_0() { return cLeftCurlyBracketKeyword_4_0; }
		
		//ports=Ports? & labels?=Labels? & size=Size? & removeIsolated?='remove isolated'?
		public UnorderedGroup getUnorderedGroup_4_1() { return cUnorderedGroup_4_1; }
		
		//ports=Ports?
		public Assignment getPortsAssignment_4_1_0() { return cPortsAssignment_4_1_0; }
		
		//Ports
		public RuleCall getPortsPortsParserRuleCall_4_1_0_0() { return cPortsPortsParserRuleCall_4_1_0_0; }
		
		//labels?=Labels?
		public Assignment getLabelsAssignment_4_1_1() { return cLabelsAssignment_4_1_1; }
		
		//Labels
		public RuleCall getLabelsLabelsParserRuleCall_4_1_1_0() { return cLabelsLabelsParserRuleCall_4_1_1_0; }
		
		//size=Size?
		public Assignment getSizeAssignment_4_1_2() { return cSizeAssignment_4_1_2; }
		
		//Size
		public RuleCall getSizeSizeParserRuleCall_4_1_2_0() { return cSizeSizeParserRuleCall_4_1_2_0; }
		
		//removeIsolated?='remove isolated'?
		public Assignment getRemoveIsolatedAssignment_4_1_3() { return cRemoveIsolatedAssignment_4_1_3; }
		
		//'remove isolated'
		public Keyword getRemoveIsolatedRemoveIsolatedKeyword_4_1_3_0() { return cRemoveIsolatedRemoveIsolatedKeyword_4_1_3_0; }
		
		//'}'
		public Keyword getRightCurlyBracketKeyword_4_2() { return cRightCurlyBracketKeyword_4_2; }
	}
	public class SizeElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.core.debug.grandom.GRandom.Size");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Action cSizeAction_0 = (Action)cGroup.eContents().get(0);
		private final Group cGroup_1 = (Group)cGroup.eContents().get(1);
		private final Keyword cSizeKeyword_1_0 = (Keyword)cGroup_1.eContents().get(0);
		private final Group cGroup_1_1 = (Group)cGroup_1.eContents().get(1);
		private final Keyword cLeftCurlyBracketKeyword_1_1_0 = (Keyword)cGroup_1_1.eContents().get(0);
		private final UnorderedGroup cUnorderedGroup_1_1_1 = (UnorderedGroup)cGroup_1_1.eContents().get(1);
		private final Group cGroup_1_1_1_0 = (Group)cUnorderedGroup_1_1_1.eContents().get(0);
		private final Keyword cHeightKeyword_1_1_1_0_0 = (Keyword)cGroup_1_1_1_0.eContents().get(0);
		private final Keyword cEqualsSignKeyword_1_1_1_0_1 = (Keyword)cGroup_1_1_1_0.eContents().get(1);
		private final Assignment cHeightAssignment_1_1_1_0_2 = (Assignment)cGroup_1_1_1_0.eContents().get(2);
		private final RuleCall cHeightDoubleQuantityParserRuleCall_1_1_1_0_2_0 = (RuleCall)cHeightAssignment_1_1_1_0_2.eContents().get(0);
		private final Group cGroup_1_1_1_1 = (Group)cUnorderedGroup_1_1_1.eContents().get(1);
		private final Keyword cWidthKeyword_1_1_1_1_0 = (Keyword)cGroup_1_1_1_1.eContents().get(0);
		private final Keyword cEqualsSignKeyword_1_1_1_1_1 = (Keyword)cGroup_1_1_1_1.eContents().get(1);
		private final Assignment cWidthAssignment_1_1_1_1_2 = (Assignment)cGroup_1_1_1_1.eContents().get(2);
		private final RuleCall cWidthDoubleQuantityParserRuleCall_1_1_1_1_2_0 = (RuleCall)cWidthAssignment_1_1_1_1_2.eContents().get(0);
		private final Keyword cRightCurlyBracketKeyword_1_1_2 = (Keyword)cGroup_1_1.eContents().get(2);
		
		//Size:
		//	{Size} ('size' ('{' (('height' '=' height=DoubleQuantity)? & ('width' '=' width=DoubleQuantity)?)
		//	'}')?);
		@Override public ParserRule getRule() { return rule; }
		
		//{Size} ('size' ('{' (('height' '=' height=DoubleQuantity)? & ('width' '=' width=DoubleQuantity)?) '}')?)
		public Group getGroup() { return cGroup; }
		
		//{Size}
		public Action getSizeAction_0() { return cSizeAction_0; }
		
		//'size' ('{' (('height' '=' height=DoubleQuantity)? & ('width' '=' width=DoubleQuantity)?) '}')?
		public Group getGroup_1() { return cGroup_1; }
		
		//'size'
		public Keyword getSizeKeyword_1_0() { return cSizeKeyword_1_0; }
		
		//('{' (('height' '=' height=DoubleQuantity)? & ('width' '=' width=DoubleQuantity)?) '}')?
		public Group getGroup_1_1() { return cGroup_1_1; }
		
		//'{'
		public Keyword getLeftCurlyBracketKeyword_1_1_0() { return cLeftCurlyBracketKeyword_1_1_0; }
		
		//('height' '=' height=DoubleQuantity)? & ('width' '=' width=DoubleQuantity)?
		public UnorderedGroup getUnorderedGroup_1_1_1() { return cUnorderedGroup_1_1_1; }
		
		//('height' '=' height=DoubleQuantity)?
		public Group getGroup_1_1_1_0() { return cGroup_1_1_1_0; }
		
		//'height'
		public Keyword getHeightKeyword_1_1_1_0_0() { return cHeightKeyword_1_1_1_0_0; }
		
		//'='
		public Keyword getEqualsSignKeyword_1_1_1_0_1() { return cEqualsSignKeyword_1_1_1_0_1; }
		
		//height=DoubleQuantity
		public Assignment getHeightAssignment_1_1_1_0_2() { return cHeightAssignment_1_1_1_0_2; }
		
		//DoubleQuantity
		public RuleCall getHeightDoubleQuantityParserRuleCall_1_1_1_0_2_0() { return cHeightDoubleQuantityParserRuleCall_1_1_1_0_2_0; }
		
		//('width' '=' width=DoubleQuantity)?
		public Group getGroup_1_1_1_1() { return cGroup_1_1_1_1; }
		
		//'width'
		public Keyword getWidthKeyword_1_1_1_1_0() { return cWidthKeyword_1_1_1_1_0; }
		
		//'='
		public Keyword getEqualsSignKeyword_1_1_1_1_1() { return cEqualsSignKeyword_1_1_1_1_1; }
		
		//width=DoubleQuantity
		public Assignment getWidthAssignment_1_1_1_1_2() { return cWidthAssignment_1_1_1_1_2; }
		
		//DoubleQuantity
		public RuleCall getWidthDoubleQuantityParserRuleCall_1_1_1_1_2_0() { return cWidthDoubleQuantityParserRuleCall_1_1_1_1_2_0; }
		
		//'}'
		public Keyword getRightCurlyBracketKeyword_1_1_2() { return cRightCurlyBracketKeyword_1_1_2; }
	}
	public class PortsElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.core.debug.grandom.GRandom.Ports");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Action cPortsAction_0 = (Action)cGroup.eContents().get(0);
		private final Keyword cPortsKeyword_1 = (Keyword)cGroup.eContents().get(1);
		private final Group cGroup_2 = (Group)cGroup.eContents().get(2);
		private final Keyword cLeftCurlyBracketKeyword_2_0 = (Keyword)cGroup_2.eContents().get(0);
		private final UnorderedGroup cUnorderedGroup_2_1 = (UnorderedGroup)cGroup_2.eContents().get(1);
		private final Assignment cLabelsAssignment_2_1_0 = (Assignment)cUnorderedGroup_2_1.eContents().get(0);
		private final RuleCall cLabelsLabelsParserRuleCall_2_1_0_0 = (RuleCall)cLabelsAssignment_2_1_0.eContents().get(0);
		private final Group cGroup_2_1_1 = (Group)cUnorderedGroup_2_1.eContents().get(1);
		private final Keyword cReUseKeyword_2_1_1_0 = (Keyword)cGroup_2_1_1.eContents().get(0);
		private final Keyword cEqualsSignKeyword_2_1_1_1 = (Keyword)cGroup_2_1_1.eContents().get(1);
		private final Assignment cReUseAssignment_2_1_1_2 = (Assignment)cGroup_2_1_1.eContents().get(2);
		private final RuleCall cReUseDoubleQuantityParserRuleCall_2_1_1_2_0 = (RuleCall)cReUseAssignment_2_1_1_2.eContents().get(0);
		private final Assignment cSizeAssignment_2_1_2 = (Assignment)cUnorderedGroup_2_1.eContents().get(2);
		private final RuleCall cSizeSizeParserRuleCall_2_1_2_0 = (RuleCall)cSizeAssignment_2_1_2.eContents().get(0);
		private final Group cGroup_2_1_3 = (Group)cUnorderedGroup_2_1.eContents().get(3);
		private final Keyword cConstraintKeyword_2_1_3_0 = (Keyword)cGroup_2_1_3.eContents().get(0);
		private final Keyword cEqualsSignKeyword_2_1_3_1 = (Keyword)cGroup_2_1_3.eContents().get(1);
		private final Assignment cConstraintAssignment_2_1_3_2 = (Assignment)cGroup_2_1_3.eContents().get(2);
		private final RuleCall cConstraintConstraintTypeEnumRuleCall_2_1_3_2_0 = (RuleCall)cConstraintAssignment_2_1_3_2.eContents().get(0);
		private final Assignment cFlowAssignment_2_1_4 = (Assignment)cUnorderedGroup_2_1.eContents().get(4);
		private final RuleCall cFlowFlowParserRuleCall_2_1_4_0 = (RuleCall)cFlowAssignment_2_1_4.eContents().get(0);
		private final Keyword cRightCurlyBracketKeyword_2_2 = (Keyword)cGroup_2.eContents().get(2);
		
		//Ports:
		//	{Ports}
		//	'ports' ('{' (labels?=Labels?
		//	& ('re-use' '=' reUse=DoubleQuantity)?
		//	& size=Size?
		//	& ('constraint' '=' constraint=ConstraintType)?
		//	& flow+=Flow*)
		//	'}')?;
		@Override public ParserRule getRule() { return rule; }
		
		//{Ports} 'ports' ('{' (labels?=Labels? & ('re-use' '=' reUse=DoubleQuantity)? & size=Size? & ('constraint' '='
		//constraint=ConstraintType)? & flow+=Flow*) '}')?
		public Group getGroup() { return cGroup; }
		
		//{Ports}
		public Action getPortsAction_0() { return cPortsAction_0; }
		
		//'ports'
		public Keyword getPortsKeyword_1() { return cPortsKeyword_1; }
		
		//('{' (labels?=Labels? & ('re-use' '=' reUse=DoubleQuantity)? & size=Size? & ('constraint' '='
		//constraint=ConstraintType)? & flow+=Flow*) '}')?
		public Group getGroup_2() { return cGroup_2; }
		
		//'{'
		public Keyword getLeftCurlyBracketKeyword_2_0() { return cLeftCurlyBracketKeyword_2_0; }
		
		//labels?=Labels? & ('re-use' '=' reUse=DoubleQuantity)? & size=Size? & ('constraint' '=' constraint=ConstraintType)? &
		//flow+=Flow*
		public UnorderedGroup getUnorderedGroup_2_1() { return cUnorderedGroup_2_1; }
		
		//labels?=Labels?
		public Assignment getLabelsAssignment_2_1_0() { return cLabelsAssignment_2_1_0; }
		
		//Labels
		public RuleCall getLabelsLabelsParserRuleCall_2_1_0_0() { return cLabelsLabelsParserRuleCall_2_1_0_0; }
		
		//('re-use' '=' reUse=DoubleQuantity)?
		public Group getGroup_2_1_1() { return cGroup_2_1_1; }
		
		//'re-use'
		public Keyword getReUseKeyword_2_1_1_0() { return cReUseKeyword_2_1_1_0; }
		
		//'='
		public Keyword getEqualsSignKeyword_2_1_1_1() { return cEqualsSignKeyword_2_1_1_1; }
		
		//reUse=DoubleQuantity
		public Assignment getReUseAssignment_2_1_1_2() { return cReUseAssignment_2_1_1_2; }
		
		//DoubleQuantity
		public RuleCall getReUseDoubleQuantityParserRuleCall_2_1_1_2_0() { return cReUseDoubleQuantityParserRuleCall_2_1_1_2_0; }
		
		//size=Size?
		public Assignment getSizeAssignment_2_1_2() { return cSizeAssignment_2_1_2; }
		
		//Size
		public RuleCall getSizeSizeParserRuleCall_2_1_2_0() { return cSizeSizeParserRuleCall_2_1_2_0; }
		
		//('constraint' '=' constraint=ConstraintType)?
		public Group getGroup_2_1_3() { return cGroup_2_1_3; }
		
		//'constraint'
		public Keyword getConstraintKeyword_2_1_3_0() { return cConstraintKeyword_2_1_3_0; }
		
		//'='
		public Keyword getEqualsSignKeyword_2_1_3_1() { return cEqualsSignKeyword_2_1_3_1; }
		
		//constraint=ConstraintType
		public Assignment getConstraintAssignment_2_1_3_2() { return cConstraintAssignment_2_1_3_2; }
		
		//ConstraintType
		public RuleCall getConstraintConstraintTypeEnumRuleCall_2_1_3_2_0() { return cConstraintConstraintTypeEnumRuleCall_2_1_3_2_0; }
		
		//flow+=Flow*
		public Assignment getFlowAssignment_2_1_4() { return cFlowAssignment_2_1_4; }
		
		//Flow
		public RuleCall getFlowFlowParserRuleCall_2_1_4_0() { return cFlowFlowParserRuleCall_2_1_4_0; }
		
		//'}'
		public Keyword getRightCurlyBracketKeyword_2_2() { return cRightCurlyBracketKeyword_2_2; }
	}
	public class FlowElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.core.debug.grandom.GRandom.Flow");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Assignment cFlowTypeAssignment_0 = (Assignment)cGroup.eContents().get(0);
		private final RuleCall cFlowTypeFlowTypeEnumRuleCall_0_0 = (RuleCall)cFlowTypeAssignment_0.eContents().get(0);
		private final Assignment cSideAssignment_1 = (Assignment)cGroup.eContents().get(1);
		private final RuleCall cSideSideEnumRuleCall_1_0 = (RuleCall)cSideAssignment_1.eContents().get(0);
		private final Keyword cEqualsSignKeyword_2 = (Keyword)cGroup.eContents().get(2);
		private final Assignment cAmountAssignment_3 = (Assignment)cGroup.eContents().get(3);
		private final RuleCall cAmountDoubleQuantityParserRuleCall_3_0 = (RuleCall)cAmountAssignment_3.eContents().get(0);
		
		//Flow:
		//	flowType=FlowType side=Side '=' amount=DoubleQuantity;
		@Override public ParserRule getRule() { return rule; }
		
		//flowType=FlowType side=Side '=' amount=DoubleQuantity
		public Group getGroup() { return cGroup; }
		
		//flowType=FlowType
		public Assignment getFlowTypeAssignment_0() { return cFlowTypeAssignment_0; }
		
		//FlowType
		public RuleCall getFlowTypeFlowTypeEnumRuleCall_0_0() { return cFlowTypeFlowTypeEnumRuleCall_0_0; }
		
		//side=Side
		public Assignment getSideAssignment_1() { return cSideAssignment_1; }
		
		//Side
		public RuleCall getSideSideEnumRuleCall_1_0() { return cSideSideEnumRuleCall_1_0; }
		
		//'='
		public Keyword getEqualsSignKeyword_2() { return cEqualsSignKeyword_2; }
		
		//amount=DoubleQuantity
		public Assignment getAmountAssignment_3() { return cAmountAssignment_3; }
		
		//DoubleQuantity
		public RuleCall getAmountDoubleQuantityParserRuleCall_3_0() { return cAmountDoubleQuantityParserRuleCall_3_0; }
	}
	public class LabelsElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.core.debug.grandom.GRandom.Labels");
		private final Keyword cLabelsKeyword = (Keyword)rule.eContents().get(1);
		
		//Labels:
		//	'labels';
		@Override public ParserRule getRule() { return rule; }
		
		//'labels'
		public Keyword getLabelsKeyword() { return cLabelsKeyword; }
	}
	public class DoubleQuantityElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.core.debug.grandom.GRandom.DoubleQuantity");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final Assignment cQuantAssignment_0 = (Assignment)cAlternatives.eContents().get(0);
		private final RuleCall cQuantDoubleParserRuleCall_0_0 = (RuleCall)cQuantAssignment_0.eContents().get(0);
		private final Group cGroup_1 = (Group)cAlternatives.eContents().get(1);
		private final Assignment cMinAssignment_1_0 = (Assignment)cGroup_1.eContents().get(0);
		private final RuleCall cMinDoubleParserRuleCall_1_0_0 = (RuleCall)cMinAssignment_1_0.eContents().get(0);
		private final Assignment cMinMaxAssignment_1_1 = (Assignment)cGroup_1.eContents().get(1);
		private final Keyword cMinMaxToKeyword_1_1_0 = (Keyword)cMinMaxAssignment_1_1.eContents().get(0);
		private final Assignment cMaxAssignment_1_2 = (Assignment)cGroup_1.eContents().get(2);
		private final RuleCall cMaxDoubleParserRuleCall_1_2_0 = (RuleCall)cMaxAssignment_1_2.eContents().get(0);
		private final Group cGroup_2 = (Group)cAlternatives.eContents().get(2);
		private final Assignment cMeanAssignment_2_0 = (Assignment)cGroup_2.eContents().get(0);
		private final RuleCall cMeanDoubleParserRuleCall_2_0_0 = (RuleCall)cMeanAssignment_2_0.eContents().get(0);
		private final Assignment cGaussianAssignment_2_1 = (Assignment)cGroup_2.eContents().get(1);
		private final RuleCall cGaussianPmParserRuleCall_2_1_0 = (RuleCall)cGaussianAssignment_2_1.eContents().get(0);
		private final Assignment cStddvAssignment_2_2 = (Assignment)cGroup_2.eContents().get(2);
		private final RuleCall cStddvDoubleParserRuleCall_2_2_0 = (RuleCall)cStddvAssignment_2_2.eContents().get(0);
		
		//DoubleQuantity:
		//	quant=Double | min=Double minMax?='to' max=Double | mean=Double gaussian?=Pm stddv=Double;
		@Override public ParserRule getRule() { return rule; }
		
		//quant=Double | min=Double minMax?='to' max=Double | mean=Double gaussian?=Pm stddv=Double
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//quant=Double
		public Assignment getQuantAssignment_0() { return cQuantAssignment_0; }
		
		//Double
		public RuleCall getQuantDoubleParserRuleCall_0_0() { return cQuantDoubleParserRuleCall_0_0; }
		
		//min=Double minMax?='to' max=Double
		public Group getGroup_1() { return cGroup_1; }
		
		//min=Double
		public Assignment getMinAssignment_1_0() { return cMinAssignment_1_0; }
		
		//Double
		public RuleCall getMinDoubleParserRuleCall_1_0_0() { return cMinDoubleParserRuleCall_1_0_0; }
		
		//minMax?='to'
		public Assignment getMinMaxAssignment_1_1() { return cMinMaxAssignment_1_1; }
		
		//'to'
		public Keyword getMinMaxToKeyword_1_1_0() { return cMinMaxToKeyword_1_1_0; }
		
		//max=Double
		public Assignment getMaxAssignment_1_2() { return cMaxAssignment_1_2; }
		
		//Double
		public RuleCall getMaxDoubleParserRuleCall_1_2_0() { return cMaxDoubleParserRuleCall_1_2_0; }
		
		//mean=Double gaussian?=Pm stddv=Double
		public Group getGroup_2() { return cGroup_2; }
		
		//mean=Double
		public Assignment getMeanAssignment_2_0() { return cMeanAssignment_2_0; }
		
		//Double
		public RuleCall getMeanDoubleParserRuleCall_2_0_0() { return cMeanDoubleParserRuleCall_2_0_0; }
		
		//gaussian?=Pm
		public Assignment getGaussianAssignment_2_1() { return cGaussianAssignment_2_1; }
		
		//Pm
		public RuleCall getGaussianPmParserRuleCall_2_1_0() { return cGaussianPmParserRuleCall_2_1_0; }
		
		//stddv=Double
		public Assignment getStddvAssignment_2_2() { return cStddvAssignment_2_2; }
		
		//Double
		public RuleCall getStddvDoubleParserRuleCall_2_2_0() { return cStddvDoubleParserRuleCall_2_2_0; }
	}
	public class PmElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.core.debug.grandom.GRandom.Pm");
		private final Keyword cPlusSignSolidusHyphenMinusKeyword = (Keyword)rule.eContents().get(1);
		
		//Pm:
		//	'+/-';
		@Override public ParserRule getRule() { return rule; }
		
		//'+/-'
		public Keyword getPlusSignSolidusHyphenMinusKeyword() { return cPlusSignSolidusHyphenMinusKeyword; }
	}
	public class DoubleElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.core.debug.grandom.GRandom.Double");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final RuleCall cINTTerminalRuleCall_0 = (RuleCall)cGroup.eContents().get(0);
		private final Group cGroup_1 = (Group)cGroup.eContents().get(1);
		private final Keyword cFullStopKeyword_1_0 = (Keyword)cGroup_1.eContents().get(0);
		private final RuleCall cINTTerminalRuleCall_1_1 = (RuleCall)cGroup_1.eContents().get(1);
		
		//Double ecore::EDoubleObject:
		//	INT ('.' INT)?;
		@Override public ParserRule getRule() { return rule; }
		
		//INT ('.' INT)?
		public Group getGroup() { return cGroup; }
		
		//INT
		public RuleCall getINTTerminalRuleCall_0() { return cINTTerminalRuleCall_0; }
		
		//('.' INT)?
		public Group getGroup_1() { return cGroup_1; }
		
		//'.'
		public Keyword getFullStopKeyword_1_0() { return cFullStopKeyword_1_0; }
		
		//INT
		public RuleCall getINTTerminalRuleCall_1_1() { return cINTTerminalRuleCall_1_1; }
	}
	public class FloatElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.core.debug.grandom.GRandom.Float");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final RuleCall cINTTerminalRuleCall_0 = (RuleCall)cGroup.eContents().get(0);
		private final Group cGroup_1 = (Group)cGroup.eContents().get(1);
		private final Keyword cFullStopKeyword_1_0 = (Keyword)cGroup_1.eContents().get(0);
		private final RuleCall cINTTerminalRuleCall_1_1 = (RuleCall)cGroup_1.eContents().get(1);
		
		//Float ecore::EFloatObject:
		//	INT ('.' INT)?;
		@Override public ParserRule getRule() { return rule; }
		
		//INT ('.' INT)?
		public Group getGroup() { return cGroup; }
		
		//INT
		public RuleCall getINTTerminalRuleCall_0() { return cINTTerminalRuleCall_0; }
		
		//('.' INT)?
		public Group getGroup_1() { return cGroup_1; }
		
		//'.'
		public Keyword getFullStopKeyword_1_0() { return cFullStopKeyword_1_0; }
		
		//INT
		public RuleCall getINTTerminalRuleCall_1_1() { return cINTTerminalRuleCall_1_1; }
	}
	public class IntegerElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.core.debug.grandom.GRandom.Integer");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final RuleCall cINTTerminalRuleCall_0 = (RuleCall)cGroup.eContents().get(0);
		private final Group cGroup_1 = (Group)cGroup.eContents().get(1);
		private final Keyword cFullStopKeyword_1_0 = (Keyword)cGroup_1.eContents().get(0);
		private final RuleCall cINTTerminalRuleCall_1_1 = (RuleCall)cGroup_1.eContents().get(1);
		
		//Integer ecore::EIntegerObject:
		//	INT ('.' INT)?;
		@Override public ParserRule getRule() { return rule; }
		
		//INT ('.' INT)?
		public Group getGroup() { return cGroup; }
		
		//INT
		public RuleCall getINTTerminalRuleCall_0() { return cINTTerminalRuleCall_0; }
		
		//('.' INT)?
		public Group getGroup_1() { return cGroup_1; }
		
		//'.'
		public Keyword getFullStopKeyword_1_0() { return cFullStopKeyword_1_0; }
		
		//INT
		public RuleCall getINTTerminalRuleCall_1_1() { return cINTTerminalRuleCall_1_1; }
	}
	
	public class FormatsElements extends AbstractEnumRuleElementFinder {
		private final EnumRule rule = (EnumRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.core.debug.grandom.GRandom.Formats");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final EnumLiteralDeclaration cElktEnumLiteralDeclaration_0 = (EnumLiteralDeclaration)cAlternatives.eContents().get(0);
		private final Keyword cElktElktKeyword_0_0 = (Keyword)cElktEnumLiteralDeclaration_0.eContents().get(0);
		private final EnumLiteralDeclaration cElkgEnumLiteralDeclaration_1 = (EnumLiteralDeclaration)cAlternatives.eContents().get(1);
		private final Keyword cElkgElkgKeyword_1_0 = (Keyword)cElkgEnumLiteralDeclaration_1.eContents().get(0);
		
		//enum Formats:
		//	elkt
		//	| elkg;
		public EnumRule getRule() { return rule; }
		
		//elkt | elkg
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//elkt
		public EnumLiteralDeclaration getElktEnumLiteralDeclaration_0() { return cElktEnumLiteralDeclaration_0; }
		
		//"elkt"
		public Keyword getElktElktKeyword_0_0() { return cElktElktKeyword_0_0; }
		
		//elkg
		public EnumLiteralDeclaration getElkgEnumLiteralDeclaration_1() { return cElkgEnumLiteralDeclaration_1; }
		
		//"elkg"
		public Keyword getElkgElkgKeyword_1_0() { return cElkgElkgKeyword_1_0; }
	}
	public class FormElements extends AbstractEnumRuleElementFinder {
		private final EnumRule rule = (EnumRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.core.debug.grandom.GRandom.Form");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final EnumLiteralDeclaration cTreesEnumLiteralDeclaration_0 = (EnumLiteralDeclaration)cAlternatives.eContents().get(0);
		private final Keyword cTreesTreesKeyword_0_0 = (Keyword)cTreesEnumLiteralDeclaration_0.eContents().get(0);
		private final EnumLiteralDeclaration cCustomEnumLiteralDeclaration_1 = (EnumLiteralDeclaration)cAlternatives.eContents().get(1);
		private final Keyword cCustomGraphsKeyword_1_0 = (Keyword)cCustomEnumLiteralDeclaration_1.eContents().get(0);
		private final EnumLiteralDeclaration cBipartiteEnumLiteralDeclaration_2 = (EnumLiteralDeclaration)cAlternatives.eContents().get(2);
		private final Keyword cBipartiteBipartiteGraphsKeyword_2_0 = (Keyword)cBipartiteEnumLiteralDeclaration_2.eContents().get(0);
		private final EnumLiteralDeclaration cBiconnectedEnumLiteralDeclaration_3 = (EnumLiteralDeclaration)cAlternatives.eContents().get(3);
		private final Keyword cBiconnectedBiconnectedGraphsKeyword_3_0 = (Keyword)cBiconnectedEnumLiteralDeclaration_3.eContents().get(0);
		private final EnumLiteralDeclaration cTriconnectedEnumLiteralDeclaration_4 = (EnumLiteralDeclaration)cAlternatives.eContents().get(4);
		private final Keyword cTriconnectedTriconnectedGraphsKeyword_4_0 = (Keyword)cTriconnectedEnumLiteralDeclaration_4.eContents().get(0);
		private final EnumLiteralDeclaration cAcyclicEnumLiteralDeclaration_5 = (EnumLiteralDeclaration)cAlternatives.eContents().get(5);
		private final Keyword cAcyclicAcyclicGraphsKeyword_5_0 = (Keyword)cAcyclicEnumLiteralDeclaration_5.eContents().get(0);
		
		//enum Form:
		//	trees
		//	| custom='graphs'
		//	| bipartite='bipartite graphs'
		//	| biconnected='biconnected graphs'
		//	| triconnected='triconnected graphs'
		//	| acyclic='acyclic graphs';
		public EnumRule getRule() { return rule; }
		
		//trees | custom='graphs' | bipartite='bipartite graphs' | biconnected='biconnected graphs' |
		//triconnected='triconnected graphs' | acyclic='acyclic graphs'
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//trees
		public EnumLiteralDeclaration getTreesEnumLiteralDeclaration_0() { return cTreesEnumLiteralDeclaration_0; }
		
		//"trees"
		public Keyword getTreesTreesKeyword_0_0() { return cTreesTreesKeyword_0_0; }
		
		//custom='graphs'
		public EnumLiteralDeclaration getCustomEnumLiteralDeclaration_1() { return cCustomEnumLiteralDeclaration_1; }
		
		//'graphs'
		public Keyword getCustomGraphsKeyword_1_0() { return cCustomGraphsKeyword_1_0; }
		
		//bipartite='bipartite graphs'
		public EnumLiteralDeclaration getBipartiteEnumLiteralDeclaration_2() { return cBipartiteEnumLiteralDeclaration_2; }
		
		//'bipartite graphs'
		public Keyword getBipartiteBipartiteGraphsKeyword_2_0() { return cBipartiteBipartiteGraphsKeyword_2_0; }
		
		//biconnected='biconnected graphs'
		public EnumLiteralDeclaration getBiconnectedEnumLiteralDeclaration_3() { return cBiconnectedEnumLiteralDeclaration_3; }
		
		//'biconnected graphs'
		public Keyword getBiconnectedBiconnectedGraphsKeyword_3_0() { return cBiconnectedBiconnectedGraphsKeyword_3_0; }
		
		//triconnected='triconnected graphs'
		public EnumLiteralDeclaration getTriconnectedEnumLiteralDeclaration_4() { return cTriconnectedEnumLiteralDeclaration_4; }
		
		//'triconnected graphs'
		public Keyword getTriconnectedTriconnectedGraphsKeyword_4_0() { return cTriconnectedTriconnectedGraphsKeyword_4_0; }
		
		//acyclic='acyclic graphs'
		public EnumLiteralDeclaration getAcyclicEnumLiteralDeclaration_5() { return cAcyclicEnumLiteralDeclaration_5; }
		
		//'acyclic graphs'
		public Keyword getAcyclicAcyclicGraphsKeyword_5_0() { return cAcyclicAcyclicGraphsKeyword_5_0; }
	}
	public class SideElements extends AbstractEnumRuleElementFinder {
		private final EnumRule rule = (EnumRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.core.debug.grandom.GRandom.Side");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final EnumLiteralDeclaration cNorthEnumLiteralDeclaration_0 = (EnumLiteralDeclaration)cAlternatives.eContents().get(0);
		private final Keyword cNorthNorthKeyword_0_0 = (Keyword)cNorthEnumLiteralDeclaration_0.eContents().get(0);
		private final EnumLiteralDeclaration cEastEnumLiteralDeclaration_1 = (EnumLiteralDeclaration)cAlternatives.eContents().get(1);
		private final Keyword cEastEastKeyword_1_0 = (Keyword)cEastEnumLiteralDeclaration_1.eContents().get(0);
		private final EnumLiteralDeclaration cSouthEnumLiteralDeclaration_2 = (EnumLiteralDeclaration)cAlternatives.eContents().get(2);
		private final Keyword cSouthSouthKeyword_2_0 = (Keyword)cSouthEnumLiteralDeclaration_2.eContents().get(0);
		private final EnumLiteralDeclaration cWestEnumLiteralDeclaration_3 = (EnumLiteralDeclaration)cAlternatives.eContents().get(3);
		private final Keyword cWestWestKeyword_3_0 = (Keyword)cWestEnumLiteralDeclaration_3.eContents().get(0);
		
		//enum Side:
		//	north
		//	| east
		//	| south
		//	| west;
		public EnumRule getRule() { return rule; }
		
		//north | east | south | west
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//north
		public EnumLiteralDeclaration getNorthEnumLiteralDeclaration_0() { return cNorthEnumLiteralDeclaration_0; }
		
		//"north"
		public Keyword getNorthNorthKeyword_0_0() { return cNorthNorthKeyword_0_0; }
		
		//east
		public EnumLiteralDeclaration getEastEnumLiteralDeclaration_1() { return cEastEnumLiteralDeclaration_1; }
		
		//"east"
		public Keyword getEastEastKeyword_1_0() { return cEastEastKeyword_1_0; }
		
		//south
		public EnumLiteralDeclaration getSouthEnumLiteralDeclaration_2() { return cSouthEnumLiteralDeclaration_2; }
		
		//"south"
		public Keyword getSouthSouthKeyword_2_0() { return cSouthSouthKeyword_2_0; }
		
		//west
		public EnumLiteralDeclaration getWestEnumLiteralDeclaration_3() { return cWestEnumLiteralDeclaration_3; }
		
		//"west"
		public Keyword getWestWestKeyword_3_0() { return cWestWestKeyword_3_0; }
	}
	public class FlowTypeElements extends AbstractEnumRuleElementFinder {
		private final EnumRule rule = (EnumRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.core.debug.grandom.GRandom.FlowType");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final EnumLiteralDeclaration cIncomingEnumLiteralDeclaration_0 = (EnumLiteralDeclaration)cAlternatives.eContents().get(0);
		private final Keyword cIncomingIncomingKeyword_0_0 = (Keyword)cIncomingEnumLiteralDeclaration_0.eContents().get(0);
		private final EnumLiteralDeclaration cOutgoingEnumLiteralDeclaration_1 = (EnumLiteralDeclaration)cAlternatives.eContents().get(1);
		private final Keyword cOutgoingOutgoingKeyword_1_0 = (Keyword)cOutgoingEnumLiteralDeclaration_1.eContents().get(0);
		
		//enum FlowType:
		//	incoming
		//	| outgoing;
		public EnumRule getRule() { return rule; }
		
		//incoming | outgoing
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//incoming
		public EnumLiteralDeclaration getIncomingEnumLiteralDeclaration_0() { return cIncomingEnumLiteralDeclaration_0; }
		
		//"incoming"
		public Keyword getIncomingIncomingKeyword_0_0() { return cIncomingIncomingKeyword_0_0; }
		
		//outgoing
		public EnumLiteralDeclaration getOutgoingEnumLiteralDeclaration_1() { return cOutgoingEnumLiteralDeclaration_1; }
		
		//"outgoing"
		public Keyword getOutgoingOutgoingKeyword_1_0() { return cOutgoingOutgoingKeyword_1_0; }
	}
	public class ConstraintTypeElements extends AbstractEnumRuleElementFinder {
		private final EnumRule rule = (EnumRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.core.debug.grandom.GRandom.ConstraintType");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final EnumLiteralDeclaration cFreeEnumLiteralDeclaration_0 = (EnumLiteralDeclaration)cAlternatives.eContents().get(0);
		private final Keyword cFreeFreeKeyword_0_0 = (Keyword)cFreeEnumLiteralDeclaration_0.eContents().get(0);
		private final EnumLiteralDeclaration cSideEnumLiteralDeclaration_1 = (EnumLiteralDeclaration)cAlternatives.eContents().get(1);
		private final Keyword cSideSideKeyword_1_0 = (Keyword)cSideEnumLiteralDeclaration_1.eContents().get(0);
		private final EnumLiteralDeclaration cPositionEnumLiteralDeclaration_2 = (EnumLiteralDeclaration)cAlternatives.eContents().get(2);
		private final Keyword cPositionPositionKeyword_2_0 = (Keyword)cPositionEnumLiteralDeclaration_2.eContents().get(0);
		private final EnumLiteralDeclaration cOrderEnumLiteralDeclaration_3 = (EnumLiteralDeclaration)cAlternatives.eContents().get(3);
		private final Keyword cOrderOrderKeyword_3_0 = (Keyword)cOrderEnumLiteralDeclaration_3.eContents().get(0);
		private final EnumLiteralDeclaration cRatioEnumLiteralDeclaration_4 = (EnumLiteralDeclaration)cAlternatives.eContents().get(4);
		private final Keyword cRatioRatioKeyword_4_0 = (Keyword)cRatioEnumLiteralDeclaration_4.eContents().get(0);
		
		//enum ConstraintType:
		//	free
		//	| side
		//	| position
		//	| order
		//	| ratio;
		public EnumRule getRule() { return rule; }
		
		//free | side | position | order | ratio
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//free
		public EnumLiteralDeclaration getFreeEnumLiteralDeclaration_0() { return cFreeEnumLiteralDeclaration_0; }
		
		//"free"
		public Keyword getFreeFreeKeyword_0_0() { return cFreeFreeKeyword_0_0; }
		
		//side
		public EnumLiteralDeclaration getSideEnumLiteralDeclaration_1() { return cSideEnumLiteralDeclaration_1; }
		
		//"side"
		public Keyword getSideSideKeyword_1_0() { return cSideSideKeyword_1_0; }
		
		//position
		public EnumLiteralDeclaration getPositionEnumLiteralDeclaration_2() { return cPositionEnumLiteralDeclaration_2; }
		
		//"position"
		public Keyword getPositionPositionKeyword_2_0() { return cPositionPositionKeyword_2_0; }
		
		//order
		public EnumLiteralDeclaration getOrderEnumLiteralDeclaration_3() { return cOrderEnumLiteralDeclaration_3; }
		
		//"order"
		public Keyword getOrderOrderKeyword_3_0() { return cOrderOrderKeyword_3_0; }
		
		//ratio
		public EnumLiteralDeclaration getRatioEnumLiteralDeclaration_4() { return cRatioEnumLiteralDeclaration_4; }
		
		//"ratio"
		public Keyword getRatioRatioKeyword_4_0() { return cRatioRatioKeyword_4_0; }
	}
	
	private final RandGraphElements pRandGraph;
	private final FormatsElements eFormats;
	private final ConfigurationElements pConfiguration;
	private final HierarchyElements pHierarchy;
	private final FormElements eForm;
	private final EdgesElements pEdges;
	private final NodesElements pNodes;
	private final SizeElements pSize;
	private final PortsElements pPorts;
	private final FlowElements pFlow;
	private final SideElements eSide;
	private final FlowTypeElements eFlowType;
	private final ConstraintTypeElements eConstraintType;
	private final LabelsElements pLabels;
	private final DoubleQuantityElements pDoubleQuantity;
	private final PmElements pPm;
	private final DoubleElements pDouble;
	private final FloatElements pFloat;
	private final IntegerElements pInteger;
	
	private final Grammar grammar;
	
	private final TerminalsGrammarAccess gaTerminals;

	@Inject
	public GRandomGrammarAccess(GrammarProvider grammarProvider,
			TerminalsGrammarAccess gaTerminals) {
		this.grammar = internalFindGrammar(grammarProvider);
		this.gaTerminals = gaTerminals;
		this.pRandGraph = new RandGraphElements();
		this.eFormats = new FormatsElements();
		this.pConfiguration = new ConfigurationElements();
		this.pHierarchy = new HierarchyElements();
		this.eForm = new FormElements();
		this.pEdges = new EdgesElements();
		this.pNodes = new NodesElements();
		this.pSize = new SizeElements();
		this.pPorts = new PortsElements();
		this.pFlow = new FlowElements();
		this.eSide = new SideElements();
		this.eFlowType = new FlowTypeElements();
		this.eConstraintType = new ConstraintTypeElements();
		this.pLabels = new LabelsElements();
		this.pDoubleQuantity = new DoubleQuantityElements();
		this.pPm = new PmElements();
		this.pDouble = new DoubleElements();
		this.pFloat = new FloatElements();
		this.pInteger = new IntegerElements();
	}
	
	protected Grammar internalFindGrammar(GrammarProvider grammarProvider) {
		Grammar grammar = grammarProvider.getGrammar(this);
		while (grammar != null) {
			if ("org.eclipse.elk.core.debug.grandom.GRandom".equals(grammar.getName())) {
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

	
	//RandGraph:
	//	configs+=Configuration*;
	public RandGraphElements getRandGraphAccess() {
		return pRandGraph;
	}
	
	public ParserRule getRandGraphRule() {
		return getRandGraphAccess().getRule();
	}
	
	//enum Formats:
	//	elkt
	//	| elkg;
	public FormatsElements getFormatsAccess() {
		return eFormats;
	}
	
	public EnumRule getFormatsRule() {
		return getFormatsAccess().getRule();
	}
	
	//Configuration:
	//	'generate' samples=INT
	//	form=Form ('{' (nodes=Nodes?
	//	& edges=Edges?
	//	& (mW?='maxWidth' '=' maxWidth=Integer)?
	//	& (mD?='maxDegree' '=' maxDegree=Integer)?
	//	& (pF?='partitionFraction' '=' fraction=DoubleQuantity)?
	//	& hierarchy=Hierarchy?
	//	& ('seed' '=' seed=Integer)?
	//	& ('format' '=' format=Formats)?
	//	& ('filename' '=' filename=STRING)?)
	//	'}')?;
	public ConfigurationElements getConfigurationAccess() {
		return pConfiguration;
	}
	
	public ParserRule getConfigurationRule() {
		return getConfigurationAccess().getRule();
	}
	
	//Hierarchy:
	//	{Hierarchy}
	//	'hierarchy' ('{' (('levels' '=' levels=DoubleQuantity)? // MISSING
	//	& ('cross-hierarchy edges' '=' edges=DoubleQuantity)?
	//	& ('compound nodes' '=' numHierarchNodes=DoubleQuantity)?
	//	& ('cross-hierarchy relative edges' '=' crossHierarchRel=DoubleQuantity)?)
	//	'}')?;
	public HierarchyElements getHierarchyAccess() {
		return pHierarchy;
	}
	
	public ParserRule getHierarchyRule() {
		return getHierarchyAccess().getRule();
	}
	
	//enum Form:
	//	trees
	//	| custom='graphs'
	//	| bipartite='bipartite graphs'
	//	| biconnected='biconnected graphs'
	//	| triconnected='triconnected graphs'
	//	| acyclic='acyclic graphs';
	public FormElements getFormAccess() {
		return eForm;
	}
	
	public EnumRule getFormRule() {
		return getFormAccess().getRule();
	}
	
	//Edges:
	//	('edges' (density?='density' | total?='total' | relative?='relative' | outbound?='outgoing')
	//	'=' nEdges=DoubleQuantity) ('{' (labels?='labels'? & selfLoops?='self loops'?)
	//	'}')?;
	public EdgesElements getEdgesAccess() {
		return pEdges;
	}
	
	public ParserRule getEdgesRule() {
		return getEdgesAccess().getRule();
	}
	
	//Nodes:
	//	{Nodes}
	//	'nodes' '=' nNodes=DoubleQuantity ('{' (ports=Ports? & labels?=Labels? & size=Size? &
	//	removeIsolated?='remove isolated'?)
	//	'}')?;
	public NodesElements getNodesAccess() {
		return pNodes;
	}
	
	public ParserRule getNodesRule() {
		return getNodesAccess().getRule();
	}
	
	//Size:
	//	{Size} ('size' ('{' (('height' '=' height=DoubleQuantity)? & ('width' '=' width=DoubleQuantity)?)
	//	'}')?);
	public SizeElements getSizeAccess() {
		return pSize;
	}
	
	public ParserRule getSizeRule() {
		return getSizeAccess().getRule();
	}
	
	//Ports:
	//	{Ports}
	//	'ports' ('{' (labels?=Labels?
	//	& ('re-use' '=' reUse=DoubleQuantity)?
	//	& size=Size?
	//	& ('constraint' '=' constraint=ConstraintType)?
	//	& flow+=Flow*)
	//	'}')?;
	public PortsElements getPortsAccess() {
		return pPorts;
	}
	
	public ParserRule getPortsRule() {
		return getPortsAccess().getRule();
	}
	
	//Flow:
	//	flowType=FlowType side=Side '=' amount=DoubleQuantity;
	public FlowElements getFlowAccess() {
		return pFlow;
	}
	
	public ParserRule getFlowRule() {
		return getFlowAccess().getRule();
	}
	
	//enum Side:
	//	north
	//	| east
	//	| south
	//	| west;
	public SideElements getSideAccess() {
		return eSide;
	}
	
	public EnumRule getSideRule() {
		return getSideAccess().getRule();
	}
	
	//enum FlowType:
	//	incoming
	//	| outgoing;
	public FlowTypeElements getFlowTypeAccess() {
		return eFlowType;
	}
	
	public EnumRule getFlowTypeRule() {
		return getFlowTypeAccess().getRule();
	}
	
	//enum ConstraintType:
	//	free
	//	| side
	//	| position
	//	| order
	//	| ratio;
	public ConstraintTypeElements getConstraintTypeAccess() {
		return eConstraintType;
	}
	
	public EnumRule getConstraintTypeRule() {
		return getConstraintTypeAccess().getRule();
	}
	
	//Labels:
	//	'labels';
	public LabelsElements getLabelsAccess() {
		return pLabels;
	}
	
	public ParserRule getLabelsRule() {
		return getLabelsAccess().getRule();
	}
	
	//DoubleQuantity:
	//	quant=Double | min=Double minMax?='to' max=Double | mean=Double gaussian?=Pm stddv=Double;
	public DoubleQuantityElements getDoubleQuantityAccess() {
		return pDoubleQuantity;
	}
	
	public ParserRule getDoubleQuantityRule() {
		return getDoubleQuantityAccess().getRule();
	}
	
	//Pm:
	//	'+/-';
	public PmElements getPmAccess() {
		return pPm;
	}
	
	public ParserRule getPmRule() {
		return getPmAccess().getRule();
	}
	
	//Double ecore::EDoubleObject:
	//	INT ('.' INT)?;
	public DoubleElements getDoubleAccess() {
		return pDouble;
	}
	
	public ParserRule getDoubleRule() {
		return getDoubleAccess().getRule();
	}
	
	//Float ecore::EFloatObject:
	//	INT ('.' INT)?;
	public FloatElements getFloatAccess() {
		return pFloat;
	}
	
	public ParserRule getFloatRule() {
		return getFloatAccess().getRule();
	}
	
	//Integer ecore::EIntegerObject:
	//	INT ('.' INT)?;
	public IntegerElements getIntegerAccess() {
		return pInteger;
	}
	
	public ParserRule getIntegerRule() {
		return getIntegerAccess().getRule();
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
