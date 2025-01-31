/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.core.meta.services;

import com.google.inject.Inject;
import com.google.inject.Singleton;
import java.util.List;
import org.eclipse.xtext.Action;
import org.eclipse.xtext.Alternatives;
import org.eclipse.xtext.Assignment;
import org.eclipse.xtext.CrossReference;
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
import org.eclipse.xtext.service.AbstractElementFinder.AbstractGrammarElementFinder;
import org.eclipse.xtext.service.GrammarProvider;
import org.eclipse.xtext.xbase.services.XbaseGrammarAccess;
import org.eclipse.xtext.xbase.services.XtypeGrammarAccess;

@Singleton
public class MetaDataGrammarAccess extends AbstractGrammarElementFinder {
	
	public class MdModelElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.core.meta.MetaData.MdModel");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Keyword cPackageKeyword_0 = (Keyword)cGroup.eContents().get(0);
		private final Assignment cNameAssignment_1 = (Assignment)cGroup.eContents().get(1);
		private final RuleCall cNameQualifiedNameParserRuleCall_1_0 = (RuleCall)cNameAssignment_1.eContents().get(0);
		private final Assignment cImportSectionAssignment_2 = (Assignment)cGroup.eContents().get(2);
		private final RuleCall cImportSectionXImportSectionParserRuleCall_2_0 = (RuleCall)cImportSectionAssignment_2.eContents().get(0);
		private final Assignment cBundleAssignment_3 = (Assignment)cGroup.eContents().get(3);
		private final RuleCall cBundleMdBundleParserRuleCall_3_0 = (RuleCall)cBundleAssignment_3.eContents().get(0);
		
		//MdModel:
		//	('package' name=QualifiedName
		//	importSection=XImportSection?
		//	bundle=MdBundle)?;
		@Override public ParserRule getRule() { return rule; }
		
		//('package' name=QualifiedName importSection=XImportSection? bundle=MdBundle)?
		public Group getGroup() { return cGroup; }
		
		//'package'
		public Keyword getPackageKeyword_0() { return cPackageKeyword_0; }
		
		//name=QualifiedName
		public Assignment getNameAssignment_1() { return cNameAssignment_1; }
		
		//QualifiedName
		public RuleCall getNameQualifiedNameParserRuleCall_1_0() { return cNameQualifiedNameParserRuleCall_1_0; }
		
		//importSection=XImportSection?
		public Assignment getImportSectionAssignment_2() { return cImportSectionAssignment_2; }
		
		//XImportSection
		public RuleCall getImportSectionXImportSectionParserRuleCall_2_0() { return cImportSectionXImportSectionParserRuleCall_2_0; }
		
		//bundle=MdBundle
		public Assignment getBundleAssignment_3() { return cBundleAssignment_3; }
		
		//MdBundle
		public RuleCall getBundleMdBundleParserRuleCall_3_0() { return cBundleMdBundleParserRuleCall_3_0; }
	}
	public class MdBundleElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.core.meta.MetaData.MdBundle");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Action cMdBundleAction_0 = (Action)cGroup.eContents().get(0);
		private final Group cGroup_1 = (Group)cGroup.eContents().get(1);
		private final Keyword cBundleKeyword_1_0 = (Keyword)cGroup_1.eContents().get(0);
		private final Keyword cLeftCurlyBracketKeyword_1_1 = (Keyword)cGroup_1.eContents().get(1);
		private final UnorderedGroup cUnorderedGroup_1_2 = (UnorderedGroup)cGroup_1.eContents().get(2);
		private final Group cGroup_1_2_0 = (Group)cUnorderedGroup_1_2.eContents().get(0);
		private final Keyword cLabelKeyword_1_2_0_0 = (Keyword)cGroup_1_2_0.eContents().get(0);
		private final Assignment cLabelAssignment_1_2_0_1 = (Assignment)cGroup_1_2_0.eContents().get(1);
		private final RuleCall cLabelSTRINGTerminalRuleCall_1_2_0_1_0 = (RuleCall)cLabelAssignment_1_2_0_1.eContents().get(0);
		private final Group cGroup_1_2_1 = (Group)cUnorderedGroup_1_2.eContents().get(1);
		private final Keyword cMetadataClassKeyword_1_2_1_0 = (Keyword)cGroup_1_2_1.eContents().get(0);
		private final Assignment cTargetClassAssignment_1_2_1_1 = (Assignment)cGroup_1_2_1.eContents().get(1);
		private final RuleCall cTargetClassQualifiedNameParserRuleCall_1_2_1_1_0 = (RuleCall)cTargetClassAssignment_1_2_1_1.eContents().get(0);
		private final Group cGroup_1_2_2 = (Group)cUnorderedGroup_1_2.eContents().get(2);
		private final Keyword cDocumentationFolderKeyword_1_2_2_0 = (Keyword)cGroup_1_2_2.eContents().get(0);
		private final Assignment cDocumentationFolderAssignment_1_2_2_1 = (Assignment)cGroup_1_2_2.eContents().get(1);
		private final RuleCall cDocumentationFolderPathParserRuleCall_1_2_2_1_0 = (RuleCall)cDocumentationFolderAssignment_1_2_2_1.eContents().get(0);
		private final Group cGroup_1_2_3 = (Group)cUnorderedGroup_1_2.eContents().get(3);
		private final Keyword cIdPrefixKeyword_1_2_3_0 = (Keyword)cGroup_1_2_3.eContents().get(0);
		private final Assignment cIdPrefixAssignment_1_2_3_1 = (Assignment)cGroup_1_2_3.eContents().get(1);
		private final RuleCall cIdPrefixQualifiedNameParserRuleCall_1_2_3_1_0 = (RuleCall)cIdPrefixAssignment_1_2_3_1.eContents().get(0);
		private final Keyword cRightCurlyBracketKeyword_1_3 = (Keyword)cGroup_1.eContents().get(3);
		private final Assignment cMembersAssignment_2 = (Assignment)cGroup.eContents().get(2);
		private final RuleCall cMembersMdBundleMemberParserRuleCall_2_0 = (RuleCall)cMembersAssignment_2.eContents().get(0);
		
		//MdBundle:
		//	{MdBundle} ('bundle' '{' (('label' label=STRING)?
		//	& ('metadataClass' targetClass=QualifiedName)?
		//	& ('documentationFolder' documentationFolder=Path)?
		//	& ('idPrefix' idPrefix=QualifiedName)?)
		//	'}')?
		//	members+=MdBundleMember*;
		@Override public ParserRule getRule() { return rule; }
		
		//{MdBundle} ('bundle' '{' (('label' label=STRING)? & ('metadataClass' targetClass=QualifiedName)? &
		//('documentationFolder' documentationFolder=Path)? & ('idPrefix' idPrefix=QualifiedName)?) '}')?
		//members+=MdBundleMember*
		public Group getGroup() { return cGroup; }
		
		//{MdBundle}
		public Action getMdBundleAction_0() { return cMdBundleAction_0; }
		
		//('bundle' '{' (('label' label=STRING)? & ('metadataClass' targetClass=QualifiedName)? & ('documentationFolder'
		//documentationFolder=Path)? & ('idPrefix' idPrefix=QualifiedName)?) '}')?
		public Group getGroup_1() { return cGroup_1; }
		
		//'bundle'
		public Keyword getBundleKeyword_1_0() { return cBundleKeyword_1_0; }
		
		//'{'
		public Keyword getLeftCurlyBracketKeyword_1_1() { return cLeftCurlyBracketKeyword_1_1; }
		
		//('label' label=STRING)? & ('metadataClass' targetClass=QualifiedName)? & ('documentationFolder'
		//documentationFolder=Path)? & ('idPrefix' idPrefix=QualifiedName)?
		public UnorderedGroup getUnorderedGroup_1_2() { return cUnorderedGroup_1_2; }
		
		//('label' label=STRING)?
		public Group getGroup_1_2_0() { return cGroup_1_2_0; }
		
		//'label'
		public Keyword getLabelKeyword_1_2_0_0() { return cLabelKeyword_1_2_0_0; }
		
		//label=STRING
		public Assignment getLabelAssignment_1_2_0_1() { return cLabelAssignment_1_2_0_1; }
		
		//STRING
		public RuleCall getLabelSTRINGTerminalRuleCall_1_2_0_1_0() { return cLabelSTRINGTerminalRuleCall_1_2_0_1_0; }
		
		//('metadataClass' targetClass=QualifiedName)?
		public Group getGroup_1_2_1() { return cGroup_1_2_1; }
		
		//'metadataClass'
		public Keyword getMetadataClassKeyword_1_2_1_0() { return cMetadataClassKeyword_1_2_1_0; }
		
		//targetClass=QualifiedName
		public Assignment getTargetClassAssignment_1_2_1_1() { return cTargetClassAssignment_1_2_1_1; }
		
		//QualifiedName
		public RuleCall getTargetClassQualifiedNameParserRuleCall_1_2_1_1_0() { return cTargetClassQualifiedNameParserRuleCall_1_2_1_1_0; }
		
		//('documentationFolder' documentationFolder=Path)?
		public Group getGroup_1_2_2() { return cGroup_1_2_2; }
		
		//'documentationFolder'
		public Keyword getDocumentationFolderKeyword_1_2_2_0() { return cDocumentationFolderKeyword_1_2_2_0; }
		
		//documentationFolder=Path
		public Assignment getDocumentationFolderAssignment_1_2_2_1() { return cDocumentationFolderAssignment_1_2_2_1; }
		
		//Path
		public RuleCall getDocumentationFolderPathParserRuleCall_1_2_2_1_0() { return cDocumentationFolderPathParserRuleCall_1_2_2_1_0; }
		
		//('idPrefix' idPrefix=QualifiedName)?
		public Group getGroup_1_2_3() { return cGroup_1_2_3; }
		
		//'idPrefix'
		public Keyword getIdPrefixKeyword_1_2_3_0() { return cIdPrefixKeyword_1_2_3_0; }
		
		//idPrefix=QualifiedName
		public Assignment getIdPrefixAssignment_1_2_3_1() { return cIdPrefixAssignment_1_2_3_1; }
		
		//QualifiedName
		public RuleCall getIdPrefixQualifiedNameParserRuleCall_1_2_3_1_0() { return cIdPrefixQualifiedNameParserRuleCall_1_2_3_1_0; }
		
		//'}'
		public Keyword getRightCurlyBracketKeyword_1_3() { return cRightCurlyBracketKeyword_1_3; }
		
		//members+=MdBundleMember*
		public Assignment getMembersAssignment_2() { return cMembersAssignment_2; }
		
		//MdBundleMember
		public RuleCall getMembersMdBundleMemberParserRuleCall_2_0() { return cMembersMdBundleMemberParserRuleCall_2_0; }
	}
	public class MdBundleMemberElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.core.meta.MetaData.MdBundleMember");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final RuleCall cMdGroupOrOptionParserRuleCall_0 = (RuleCall)cAlternatives.eContents().get(0);
		private final RuleCall cMdAlgorithmParserRuleCall_1 = (RuleCall)cAlternatives.eContents().get(1);
		private final RuleCall cMdCategoryParserRuleCall_2 = (RuleCall)cAlternatives.eContents().get(2);
		
		//MdBundleMember:
		//	MdGroupOrOption | MdAlgorithm | MdCategory;
		@Override public ParserRule getRule() { return rule; }
		
		//MdGroupOrOption | MdAlgorithm | MdCategory
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//MdGroupOrOption
		public RuleCall getMdGroupOrOptionParserRuleCall_0() { return cMdGroupOrOptionParserRuleCall_0; }
		
		//MdAlgorithm
		public RuleCall getMdAlgorithmParserRuleCall_1() { return cMdAlgorithmParserRuleCall_1; }
		
		//MdCategory
		public RuleCall getMdCategoryParserRuleCall_2() { return cMdCategoryParserRuleCall_2; }
	}
	public class MdGroupOrOptionElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.core.meta.MetaData.MdGroupOrOption");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final RuleCall cMdGroupParserRuleCall_0 = (RuleCall)cAlternatives.eContents().get(0);
		private final RuleCall cMdOptionParserRuleCall_1 = (RuleCall)cAlternatives.eContents().get(1);
		
		//MdGroupOrOption:
		//	MdGroup | MdOption;
		@Override public ParserRule getRule() { return rule; }
		
		//MdGroup | MdOption
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//MdGroup
		public RuleCall getMdGroupParserRuleCall_0() { return cMdGroupParserRuleCall_0; }
		
		//MdOption
		public RuleCall getMdOptionParserRuleCall_1() { return cMdOptionParserRuleCall_1; }
	}
	public class MdGroupElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.core.meta.MetaData.MdGroup");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Keyword cGroupKeyword_0 = (Keyword)cGroup.eContents().get(0);
		private final Assignment cNameAssignment_1 = (Assignment)cGroup.eContents().get(1);
		private final RuleCall cNameIDTerminalRuleCall_1_0 = (RuleCall)cNameAssignment_1.eContents().get(0);
		private final Keyword cLeftCurlyBracketKeyword_2 = (Keyword)cGroup.eContents().get(2);
		private final Group cGroup_3 = (Group)cGroup.eContents().get(3);
		private final Keyword cDocumentationKeyword_3_0 = (Keyword)cGroup_3.eContents().get(0);
		private final Assignment cDocumentationAssignment_3_1 = (Assignment)cGroup_3.eContents().get(1);
		private final RuleCall cDocumentationSTRINGTerminalRuleCall_3_1_0 = (RuleCall)cDocumentationAssignment_3_1.eContents().get(0);
		private final Assignment cChildrenAssignment_4 = (Assignment)cGroup.eContents().get(4);
		private final RuleCall cChildrenMdGroupOrOptionParserRuleCall_4_0 = (RuleCall)cChildrenAssignment_4.eContents().get(0);
		private final Keyword cRightCurlyBracketKeyword_5 = (Keyword)cGroup.eContents().get(5);
		
		//MdGroup:
		//	'group' name=ID '{' ('documentation' documentation=STRING)?
		//	children+=MdGroupOrOption*
		//	'}';
		@Override public ParserRule getRule() { return rule; }
		
		//'group' name=ID '{' ('documentation' documentation=STRING)? children+=MdGroupOrOption* '}'
		public Group getGroup() { return cGroup; }
		
		//'group'
		public Keyword getGroupKeyword_0() { return cGroupKeyword_0; }
		
		//name=ID
		public Assignment getNameAssignment_1() { return cNameAssignment_1; }
		
		//ID
		public RuleCall getNameIDTerminalRuleCall_1_0() { return cNameIDTerminalRuleCall_1_0; }
		
		//'{'
		public Keyword getLeftCurlyBracketKeyword_2() { return cLeftCurlyBracketKeyword_2; }
		
		//('documentation' documentation=STRING)?
		public Group getGroup_3() { return cGroup_3; }
		
		//'documentation'
		public Keyword getDocumentationKeyword_3_0() { return cDocumentationKeyword_3_0; }
		
		//documentation=STRING
		public Assignment getDocumentationAssignment_3_1() { return cDocumentationAssignment_3_1; }
		
		//STRING
		public RuleCall getDocumentationSTRINGTerminalRuleCall_3_1_0() { return cDocumentationSTRINGTerminalRuleCall_3_1_0; }
		
		//children+=MdGroupOrOption*
		public Assignment getChildrenAssignment_4() { return cChildrenAssignment_4; }
		
		//MdGroupOrOption
		public RuleCall getChildrenMdGroupOrOptionParserRuleCall_4_0() { return cChildrenMdGroupOrOptionParserRuleCall_4_0; }
		
		//'}'
		public Keyword getRightCurlyBracketKeyword_5() { return cRightCurlyBracketKeyword_5; }
	}
	public class MdOptionElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.core.meta.MetaData.MdOption");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Assignment cDeprecatedAssignment_0 = (Assignment)cGroup.eContents().get(0);
		private final Keyword cDeprecatedDeprecatedKeyword_0_0 = (Keyword)cDeprecatedAssignment_0.eContents().get(0);
		private final Alternatives cAlternatives_1 = (Alternatives)cGroup.eContents().get(1);
		private final Assignment cAdvancedAssignment_1_0 = (Assignment)cAlternatives_1.eContents().get(0);
		private final Keyword cAdvancedAdvancedKeyword_1_0_0 = (Keyword)cAdvancedAssignment_1_0.eContents().get(0);
		private final Assignment cProgrammaticAssignment_1_1 = (Assignment)cAlternatives_1.eContents().get(1);
		private final Keyword cProgrammaticProgrammaticKeyword_1_1_0 = (Keyword)cProgrammaticAssignment_1_1.eContents().get(0);
		private final Assignment cOutputAssignment_1_2 = (Assignment)cAlternatives_1.eContents().get(2);
		private final Keyword cOutputOutputKeyword_1_2_0 = (Keyword)cOutputAssignment_1_2.eContents().get(0);
		private final Assignment cGlobalAssignment_1_3 = (Assignment)cAlternatives_1.eContents().get(3);
		private final Keyword cGlobalGlobalKeyword_1_3_0 = (Keyword)cGlobalAssignment_1_3.eContents().get(0);
		private final Keyword cOptionKeyword_2 = (Keyword)cGroup.eContents().get(2);
		private final Assignment cNameAssignment_3 = (Assignment)cGroup.eContents().get(3);
		private final RuleCall cNameIDTerminalRuleCall_3_0 = (RuleCall)cNameAssignment_3.eContents().get(0);
		private final Group cGroup_4 = (Group)cGroup.eContents().get(4);
		private final Keyword cColonKeyword_4_0 = (Keyword)cGroup_4.eContents().get(0);
		private final Assignment cTypeAssignment_4_1 = (Assignment)cGroup_4.eContents().get(1);
		private final RuleCall cTypeJvmTypeReferenceParserRuleCall_4_1_0 = (RuleCall)cTypeAssignment_4_1.eContents().get(0);
		private final Keyword cLeftCurlyBracketKeyword_5 = (Keyword)cGroup.eContents().get(5);
		private final UnorderedGroup cUnorderedGroup_6 = (UnorderedGroup)cGroup.eContents().get(6);
		private final Group cGroup_6_0 = (Group)cUnorderedGroup_6.eContents().get(0);
		private final Keyword cLabelKeyword_6_0_0 = (Keyword)cGroup_6_0.eContents().get(0);
		private final Assignment cLabelAssignment_6_0_1 = (Assignment)cGroup_6_0.eContents().get(1);
		private final RuleCall cLabelSTRINGTerminalRuleCall_6_0_1_0 = (RuleCall)cLabelAssignment_6_0_1.eContents().get(0);
		private final Group cGroup_6_1 = (Group)cUnorderedGroup_6.eContents().get(1);
		private final Keyword cDescriptionKeyword_6_1_0 = (Keyword)cGroup_6_1.eContents().get(0);
		private final Assignment cDescriptionAssignment_6_1_1 = (Assignment)cGroup_6_1.eContents().get(1);
		private final RuleCall cDescriptionSTRINGTerminalRuleCall_6_1_1_0 = (RuleCall)cDescriptionAssignment_6_1_1.eContents().get(0);
		private final Group cGroup_6_2 = (Group)cUnorderedGroup_6.eContents().get(2);
		private final Keyword cDocumentationKeyword_6_2_0 = (Keyword)cGroup_6_2.eContents().get(0);
		private final Assignment cDocumentationAssignment_6_2_1 = (Assignment)cGroup_6_2.eContents().get(1);
		private final RuleCall cDocumentationSTRINGTerminalRuleCall_6_2_1_0 = (RuleCall)cDocumentationAssignment_6_2_1.eContents().get(0);
		private final Group cGroup_6_3 = (Group)cUnorderedGroup_6.eContents().get(3);
		private final Keyword cDefaultKeyword_6_3_0 = (Keyword)cGroup_6_3.eContents().get(0);
		private final Keyword cEqualsSignKeyword_6_3_1 = (Keyword)cGroup_6_3.eContents().get(1);
		private final Assignment cDefaultValueAssignment_6_3_2 = (Assignment)cGroup_6_3.eContents().get(2);
		private final RuleCall cDefaultValueXExpressionParserRuleCall_6_3_2_0 = (RuleCall)cDefaultValueAssignment_6_3_2.eContents().get(0);
		private final Group cGroup_6_4 = (Group)cUnorderedGroup_6.eContents().get(4);
		private final Keyword cLowerBoundKeyword_6_4_0 = (Keyword)cGroup_6_4.eContents().get(0);
		private final Keyword cEqualsSignKeyword_6_4_1 = (Keyword)cGroup_6_4.eContents().get(1);
		private final Assignment cLowerBoundAssignment_6_4_2 = (Assignment)cGroup_6_4.eContents().get(2);
		private final RuleCall cLowerBoundXExpressionParserRuleCall_6_4_2_0 = (RuleCall)cLowerBoundAssignment_6_4_2.eContents().get(0);
		private final Group cGroup_6_5 = (Group)cUnorderedGroup_6.eContents().get(5);
		private final Keyword cUpperBoundKeyword_6_5_0 = (Keyword)cGroup_6_5.eContents().get(0);
		private final Keyword cEqualsSignKeyword_6_5_1 = (Keyword)cGroup_6_5.eContents().get(1);
		private final Assignment cUpperBoundAssignment_6_5_2 = (Assignment)cGroup_6_5.eContents().get(2);
		private final RuleCall cUpperBoundXExpressionParserRuleCall_6_5_2_0 = (RuleCall)cUpperBoundAssignment_6_5_2.eContents().get(0);
		private final Group cGroup_6_6 = (Group)cUnorderedGroup_6.eContents().get(6);
		private final Keyword cTargetsKeyword_6_6_0 = (Keyword)cGroup_6_6.eContents().get(0);
		private final Assignment cTargetsAssignment_6_6_1 = (Assignment)cGroup_6_6.eContents().get(1);
		private final RuleCall cTargetsMdOptionTargetTypeEnumRuleCall_6_6_1_0 = (RuleCall)cTargetsAssignment_6_6_1.eContents().get(0);
		private final Group cGroup_6_6_2 = (Group)cGroup_6_6.eContents().get(2);
		private final Keyword cCommaKeyword_6_6_2_0 = (Keyword)cGroup_6_6_2.eContents().get(0);
		private final Assignment cTargetsAssignment_6_6_2_1 = (Assignment)cGroup_6_6_2.eContents().get(1);
		private final RuleCall cTargetsMdOptionTargetTypeEnumRuleCall_6_6_2_1_0 = (RuleCall)cTargetsAssignment_6_6_2_1.eContents().get(0);
		private final Group cGroup_6_7 = (Group)cUnorderedGroup_6.eContents().get(7);
		private final Keyword cLegacyIdsKeyword_6_7_0 = (Keyword)cGroup_6_7.eContents().get(0);
		private final Assignment cLegacyIdsAssignment_6_7_1 = (Assignment)cGroup_6_7.eContents().get(1);
		private final RuleCall cLegacyIdsQualifiedNameParserRuleCall_6_7_1_0 = (RuleCall)cLegacyIdsAssignment_6_7_1.eContents().get(0);
		private final Group cGroup_6_7_2 = (Group)cGroup_6_7.eContents().get(2);
		private final Keyword cCommaKeyword_6_7_2_0 = (Keyword)cGroup_6_7_2.eContents().get(0);
		private final Assignment cLegacyIdsAssignment_6_7_2_1 = (Assignment)cGroup_6_7_2.eContents().get(1);
		private final RuleCall cLegacyIdsQualifiedNameParserRuleCall_6_7_2_1_0 = (RuleCall)cLegacyIdsAssignment_6_7_2_1.eContents().get(0);
		private final Assignment cDependenciesAssignment_7 = (Assignment)cGroup.eContents().get(7);
		private final RuleCall cDependenciesMdOptionDependencyParserRuleCall_7_0 = (RuleCall)cDependenciesAssignment_7.eContents().get(0);
		private final Keyword cRightCurlyBracketKeyword_8 = (Keyword)cGroup.eContents().get(8);
		
		//MdOption:
		//	deprecated?='deprecated'? (advanced?='advanced' | programmatic?='programmatic' | output?='output' | global?='global')?
		//	'option' name=ID (':' type=JvmTypeReference)? '{' (('label' label=STRING)?
		//	& ('description' description=STRING)?
		//	& ('documentation' documentation=STRING)? // "@docu/priority.md"
		//	& ('default' '=' defaultValue=XExpression)?
		//	& ('lowerBound' '=' lowerBound=XExpression)?
		//	& ('upperBound' '=' upperBound=XExpression)?
		//	& ('targets' targets+=MdOptionTargetType (',' targets+=MdOptionTargetType)*)?
		//	& ('legacyIds' legacyIds+=QualifiedName (',' legacyIds+=QualifiedName)*)?) dependencies+=MdOptionDependency*
		//	'}';
		@Override public ParserRule getRule() { return rule; }
		
		//deprecated?='deprecated'? (advanced?='advanced' | programmatic?='programmatic' | output?='output' | global?='global')?
		//'option' name=ID (':' type=JvmTypeReference)? '{' (('label' label=STRING)? & ('description' description=STRING)? &
		//('documentation' documentation=STRING)? // "@docu/priority.md"
		//& ('default' '=' defaultValue=XExpression)? & ('lowerBound' '=' lowerBound=XExpression)? & ('upperBound' '='
		//upperBound=XExpression)? & ('targets' targets+=MdOptionTargetType (',' targets+=MdOptionTargetType)*)? & ('legacyIds'
		//legacyIds+=QualifiedName (',' legacyIds+=QualifiedName)*)?) dependencies+=MdOptionDependency* '}'
		public Group getGroup() { return cGroup; }
		
		//deprecated?='deprecated'?
		public Assignment getDeprecatedAssignment_0() { return cDeprecatedAssignment_0; }
		
		//'deprecated'
		public Keyword getDeprecatedDeprecatedKeyword_0_0() { return cDeprecatedDeprecatedKeyword_0_0; }
		
		//(advanced?='advanced' | programmatic?='programmatic' | output?='output' | global?='global')?
		public Alternatives getAlternatives_1() { return cAlternatives_1; }
		
		//advanced?='advanced'
		public Assignment getAdvancedAssignment_1_0() { return cAdvancedAssignment_1_0; }
		
		//'advanced'
		public Keyword getAdvancedAdvancedKeyword_1_0_0() { return cAdvancedAdvancedKeyword_1_0_0; }
		
		//programmatic?='programmatic'
		public Assignment getProgrammaticAssignment_1_1() { return cProgrammaticAssignment_1_1; }
		
		//'programmatic'
		public Keyword getProgrammaticProgrammaticKeyword_1_1_0() { return cProgrammaticProgrammaticKeyword_1_1_0; }
		
		//output?='output'
		public Assignment getOutputAssignment_1_2() { return cOutputAssignment_1_2; }
		
		//'output'
		public Keyword getOutputOutputKeyword_1_2_0() { return cOutputOutputKeyword_1_2_0; }
		
		//global?='global'
		public Assignment getGlobalAssignment_1_3() { return cGlobalAssignment_1_3; }
		
		//'global'
		public Keyword getGlobalGlobalKeyword_1_3_0() { return cGlobalGlobalKeyword_1_3_0; }
		
		//'option'
		public Keyword getOptionKeyword_2() { return cOptionKeyword_2; }
		
		//name=ID
		public Assignment getNameAssignment_3() { return cNameAssignment_3; }
		
		//ID
		public RuleCall getNameIDTerminalRuleCall_3_0() { return cNameIDTerminalRuleCall_3_0; }
		
		//(':' type=JvmTypeReference)?
		public Group getGroup_4() { return cGroup_4; }
		
		//':'
		public Keyword getColonKeyword_4_0() { return cColonKeyword_4_0; }
		
		//type=JvmTypeReference
		public Assignment getTypeAssignment_4_1() { return cTypeAssignment_4_1; }
		
		//JvmTypeReference
		public RuleCall getTypeJvmTypeReferenceParserRuleCall_4_1_0() { return cTypeJvmTypeReferenceParserRuleCall_4_1_0; }
		
		//'{'
		public Keyword getLeftCurlyBracketKeyword_5() { return cLeftCurlyBracketKeyword_5; }
		
		//('label' label=STRING)? & ('description' description=STRING)? & ('documentation' documentation=STRING)? // "@docu/priority.md"
		//& ('default' '=' defaultValue=XExpression)? & ('lowerBound' '=' lowerBound=XExpression)? & ('upperBound' '='
		//upperBound=XExpression)? & ('targets' targets+=MdOptionTargetType (',' targets+=MdOptionTargetType)*)? & ('legacyIds'
		//legacyIds+=QualifiedName (',' legacyIds+=QualifiedName)*)?
		public UnorderedGroup getUnorderedGroup_6() { return cUnorderedGroup_6; }
		
		//('label' label=STRING)?
		public Group getGroup_6_0() { return cGroup_6_0; }
		
		//'label'
		public Keyword getLabelKeyword_6_0_0() { return cLabelKeyword_6_0_0; }
		
		//label=STRING
		public Assignment getLabelAssignment_6_0_1() { return cLabelAssignment_6_0_1; }
		
		//STRING
		public RuleCall getLabelSTRINGTerminalRuleCall_6_0_1_0() { return cLabelSTRINGTerminalRuleCall_6_0_1_0; }
		
		//('description' description=STRING)?
		public Group getGroup_6_1() { return cGroup_6_1; }
		
		//'description'
		public Keyword getDescriptionKeyword_6_1_0() { return cDescriptionKeyword_6_1_0; }
		
		//description=STRING
		public Assignment getDescriptionAssignment_6_1_1() { return cDescriptionAssignment_6_1_1; }
		
		//STRING
		public RuleCall getDescriptionSTRINGTerminalRuleCall_6_1_1_0() { return cDescriptionSTRINGTerminalRuleCall_6_1_1_0; }
		
		//('documentation' documentation=STRING)?
		public Group getGroup_6_2() { return cGroup_6_2; }
		
		//'documentation'
		public Keyword getDocumentationKeyword_6_2_0() { return cDocumentationKeyword_6_2_0; }
		
		//documentation=STRING
		public Assignment getDocumentationAssignment_6_2_1() { return cDocumentationAssignment_6_2_1; }
		
		//STRING
		public RuleCall getDocumentationSTRINGTerminalRuleCall_6_2_1_0() { return cDocumentationSTRINGTerminalRuleCall_6_2_1_0; }
		
		//('default' '=' defaultValue=XExpression)?
		public Group getGroup_6_3() { return cGroup_6_3; }
		
		//'default'
		public Keyword getDefaultKeyword_6_3_0() { return cDefaultKeyword_6_3_0; }
		
		//'='
		public Keyword getEqualsSignKeyword_6_3_1() { return cEqualsSignKeyword_6_3_1; }
		
		//defaultValue=XExpression
		public Assignment getDefaultValueAssignment_6_3_2() { return cDefaultValueAssignment_6_3_2; }
		
		//XExpression
		public RuleCall getDefaultValueXExpressionParserRuleCall_6_3_2_0() { return cDefaultValueXExpressionParserRuleCall_6_3_2_0; }
		
		//('lowerBound' '=' lowerBound=XExpression)?
		public Group getGroup_6_4() { return cGroup_6_4; }
		
		//'lowerBound'
		public Keyword getLowerBoundKeyword_6_4_0() { return cLowerBoundKeyword_6_4_0; }
		
		//'='
		public Keyword getEqualsSignKeyword_6_4_1() { return cEqualsSignKeyword_6_4_1; }
		
		//lowerBound=XExpression
		public Assignment getLowerBoundAssignment_6_4_2() { return cLowerBoundAssignment_6_4_2; }
		
		//XExpression
		public RuleCall getLowerBoundXExpressionParserRuleCall_6_4_2_0() { return cLowerBoundXExpressionParserRuleCall_6_4_2_0; }
		
		//('upperBound' '=' upperBound=XExpression)?
		public Group getGroup_6_5() { return cGroup_6_5; }
		
		//'upperBound'
		public Keyword getUpperBoundKeyword_6_5_0() { return cUpperBoundKeyword_6_5_0; }
		
		//'='
		public Keyword getEqualsSignKeyword_6_5_1() { return cEqualsSignKeyword_6_5_1; }
		
		//upperBound=XExpression
		public Assignment getUpperBoundAssignment_6_5_2() { return cUpperBoundAssignment_6_5_2; }
		
		//XExpression
		public RuleCall getUpperBoundXExpressionParserRuleCall_6_5_2_0() { return cUpperBoundXExpressionParserRuleCall_6_5_2_0; }
		
		//('targets' targets+=MdOptionTargetType (',' targets+=MdOptionTargetType)*)?
		public Group getGroup_6_6() { return cGroup_6_6; }
		
		//'targets'
		public Keyword getTargetsKeyword_6_6_0() { return cTargetsKeyword_6_6_0; }
		
		//targets+=MdOptionTargetType
		public Assignment getTargetsAssignment_6_6_1() { return cTargetsAssignment_6_6_1; }
		
		//MdOptionTargetType
		public RuleCall getTargetsMdOptionTargetTypeEnumRuleCall_6_6_1_0() { return cTargetsMdOptionTargetTypeEnumRuleCall_6_6_1_0; }
		
		//(',' targets+=MdOptionTargetType)*
		public Group getGroup_6_6_2() { return cGroup_6_6_2; }
		
		//','
		public Keyword getCommaKeyword_6_6_2_0() { return cCommaKeyword_6_6_2_0; }
		
		//targets+=MdOptionTargetType
		public Assignment getTargetsAssignment_6_6_2_1() { return cTargetsAssignment_6_6_2_1; }
		
		//MdOptionTargetType
		public RuleCall getTargetsMdOptionTargetTypeEnumRuleCall_6_6_2_1_0() { return cTargetsMdOptionTargetTypeEnumRuleCall_6_6_2_1_0; }
		
		//('legacyIds' legacyIds+=QualifiedName (',' legacyIds+=QualifiedName)*)?
		public Group getGroup_6_7() { return cGroup_6_7; }
		
		//'legacyIds'
		public Keyword getLegacyIdsKeyword_6_7_0() { return cLegacyIdsKeyword_6_7_0; }
		
		//legacyIds+=QualifiedName
		public Assignment getLegacyIdsAssignment_6_7_1() { return cLegacyIdsAssignment_6_7_1; }
		
		//QualifiedName
		public RuleCall getLegacyIdsQualifiedNameParserRuleCall_6_7_1_0() { return cLegacyIdsQualifiedNameParserRuleCall_6_7_1_0; }
		
		//(',' legacyIds+=QualifiedName)*
		public Group getGroup_6_7_2() { return cGroup_6_7_2; }
		
		//','
		public Keyword getCommaKeyword_6_7_2_0() { return cCommaKeyword_6_7_2_0; }
		
		//legacyIds+=QualifiedName
		public Assignment getLegacyIdsAssignment_6_7_2_1() { return cLegacyIdsAssignment_6_7_2_1; }
		
		//QualifiedName
		public RuleCall getLegacyIdsQualifiedNameParserRuleCall_6_7_2_1_0() { return cLegacyIdsQualifiedNameParserRuleCall_6_7_2_1_0; }
		
		//dependencies+=MdOptionDependency*
		public Assignment getDependenciesAssignment_7() { return cDependenciesAssignment_7; }
		
		//MdOptionDependency
		public RuleCall getDependenciesMdOptionDependencyParserRuleCall_7_0() { return cDependenciesMdOptionDependencyParserRuleCall_7_0; }
		
		//'}'
		public Keyword getRightCurlyBracketKeyword_8() { return cRightCurlyBracketKeyword_8; }
	}
	public class MdOptionDependencyElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.core.meta.MetaData.MdOptionDependency");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Keyword cRequiresKeyword_0 = (Keyword)cGroup.eContents().get(0);
		private final Assignment cTargetAssignment_1 = (Assignment)cGroup.eContents().get(1);
		private final CrossReference cTargetMdOptionCrossReference_1_0 = (CrossReference)cTargetAssignment_1.eContents().get(0);
		private final RuleCall cTargetMdOptionQualifiedNameParserRuleCall_1_0_1 = (RuleCall)cTargetMdOptionCrossReference_1_0.eContents().get(1);
		private final Group cGroup_2 = (Group)cGroup.eContents().get(2);
		private final Keyword cEqualsSignEqualsSignKeyword_2_0 = (Keyword)cGroup_2.eContents().get(0);
		private final Assignment cValueAssignment_2_1 = (Assignment)cGroup_2.eContents().get(1);
		private final RuleCall cValueXExpressionParserRuleCall_2_1_0 = (RuleCall)cValueAssignment_2_1.eContents().get(0);
		
		//MdOptionDependency:
		//	'requires' target=[MdOption|QualifiedName] ('==' value=XExpression)?;
		@Override public ParserRule getRule() { return rule; }
		
		//'requires' target=[MdOption|QualifiedName] ('==' value=XExpression)?
		public Group getGroup() { return cGroup; }
		
		//'requires'
		public Keyword getRequiresKeyword_0() { return cRequiresKeyword_0; }
		
		//target=[MdOption|QualifiedName]
		public Assignment getTargetAssignment_1() { return cTargetAssignment_1; }
		
		//[MdOption|QualifiedName]
		public CrossReference getTargetMdOptionCrossReference_1_0() { return cTargetMdOptionCrossReference_1_0; }
		
		//QualifiedName
		public RuleCall getTargetMdOptionQualifiedNameParserRuleCall_1_0_1() { return cTargetMdOptionQualifiedNameParserRuleCall_1_0_1; }
		
		//('==' value=XExpression)?
		public Group getGroup_2() { return cGroup_2; }
		
		//'=='
		public Keyword getEqualsSignEqualsSignKeyword_2_0() { return cEqualsSignEqualsSignKeyword_2_0; }
		
		//value=XExpression
		public Assignment getValueAssignment_2_1() { return cValueAssignment_2_1; }
		
		//XExpression
		public RuleCall getValueXExpressionParserRuleCall_2_1_0() { return cValueXExpressionParserRuleCall_2_1_0; }
	}
	public class MdAlgorithmElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.core.meta.MetaData.MdAlgorithm");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Assignment cDeprecatedAssignment_0 = (Assignment)cGroup.eContents().get(0);
		private final Keyword cDeprecatedDeprecatedKeyword_0_0 = (Keyword)cDeprecatedAssignment_0.eContents().get(0);
		private final Keyword cAlgorithmKeyword_1 = (Keyword)cGroup.eContents().get(1);
		private final Assignment cNameAssignment_2 = (Assignment)cGroup.eContents().get(2);
		private final RuleCall cNameIDTerminalRuleCall_2_0 = (RuleCall)cNameAssignment_2.eContents().get(0);
		private final Keyword cLeftParenthesisKeyword_3 = (Keyword)cGroup.eContents().get(3);
		private final Assignment cProviderAssignment_4 = (Assignment)cGroup.eContents().get(4);
		private final RuleCall cProviderJvmTypeReferenceParserRuleCall_4_0 = (RuleCall)cProviderAssignment_4.eContents().get(0);
		private final Group cGroup_5 = (Group)cGroup.eContents().get(5);
		private final Keyword cNumberSignKeyword_5_0 = (Keyword)cGroup_5.eContents().get(0);
		private final Assignment cParameterAssignment_5_1 = (Assignment)cGroup_5.eContents().get(1);
		private final RuleCall cParameterIDTerminalRuleCall_5_1_0 = (RuleCall)cParameterAssignment_5_1.eContents().get(0);
		private final Keyword cRightParenthesisKeyword_6 = (Keyword)cGroup.eContents().get(6);
		private final Keyword cLeftCurlyBracketKeyword_7 = (Keyword)cGroup.eContents().get(7);
		private final UnorderedGroup cUnorderedGroup_8 = (UnorderedGroup)cGroup.eContents().get(8);
		private final Group cGroup_8_0 = (Group)cUnorderedGroup_8.eContents().get(0);
		private final Keyword cLabelKeyword_8_0_0 = (Keyword)cGroup_8_0.eContents().get(0);
		private final Assignment cLabelAssignment_8_0_1 = (Assignment)cGroup_8_0.eContents().get(1);
		private final RuleCall cLabelSTRINGTerminalRuleCall_8_0_1_0 = (RuleCall)cLabelAssignment_8_0_1.eContents().get(0);
		private final Group cGroup_8_1 = (Group)cUnorderedGroup_8.eContents().get(1);
		private final Keyword cMetadataClassKeyword_8_1_0 = (Keyword)cGroup_8_1.eContents().get(0);
		private final Assignment cTargetClassAssignment_8_1_1 = (Assignment)cGroup_8_1.eContents().get(1);
		private final RuleCall cTargetClassQualifiedNameParserRuleCall_8_1_1_0 = (RuleCall)cTargetClassAssignment_8_1_1.eContents().get(0);
		private final Group cGroup_8_2 = (Group)cUnorderedGroup_8.eContents().get(2);
		private final Keyword cDescriptionKeyword_8_2_0 = (Keyword)cGroup_8_2.eContents().get(0);
		private final Assignment cDescriptionAssignment_8_2_1 = (Assignment)cGroup_8_2.eContents().get(1);
		private final RuleCall cDescriptionSTRINGTerminalRuleCall_8_2_1_0 = (RuleCall)cDescriptionAssignment_8_2_1.eContents().get(0);
		private final Group cGroup_8_3 = (Group)cUnorderedGroup_8.eContents().get(3);
		private final Keyword cDocumentationKeyword_8_3_0 = (Keyword)cGroup_8_3.eContents().get(0);
		private final Assignment cDocumentationAssignment_8_3_1 = (Assignment)cGroup_8_3.eContents().get(1);
		private final RuleCall cDocumentationSTRINGTerminalRuleCall_8_3_1_0 = (RuleCall)cDocumentationAssignment_8_3_1.eContents().get(0);
		private final Group cGroup_8_4 = (Group)cUnorderedGroup_8.eContents().get(4);
		private final Keyword cCategoryKeyword_8_4_0 = (Keyword)cGroup_8_4.eContents().get(0);
		private final Assignment cCategoryAssignment_8_4_1 = (Assignment)cGroup_8_4.eContents().get(1);
		private final CrossReference cCategoryMdCategoryCrossReference_8_4_1_0 = (CrossReference)cCategoryAssignment_8_4_1.eContents().get(0);
		private final RuleCall cCategoryMdCategoryQualifiedNameParserRuleCall_8_4_1_0_1 = (RuleCall)cCategoryMdCategoryCrossReference_8_4_1_0.eContents().get(1);
		private final Group cGroup_8_5 = (Group)cUnorderedGroup_8.eContents().get(5);
		private final Keyword cPreviewKeyword_8_5_0 = (Keyword)cGroup_8_5.eContents().get(0);
		private final Assignment cPreviewImageAssignment_8_5_1 = (Assignment)cGroup_8_5.eContents().get(1);
		private final RuleCall cPreviewImagePathParserRuleCall_8_5_1_0 = (RuleCall)cPreviewImageAssignment_8_5_1.eContents().get(0);
		private final Group cGroup_8_6 = (Group)cUnorderedGroup_8.eContents().get(6);
		private final Keyword cFeaturesKeyword_8_6_0 = (Keyword)cGroup_8_6.eContents().get(0);
		private final Assignment cSupportedFeaturesAssignment_8_6_1 = (Assignment)cGroup_8_6.eContents().get(1);
		private final RuleCall cSupportedFeaturesMdGraphFeatureEnumRuleCall_8_6_1_0 = (RuleCall)cSupportedFeaturesAssignment_8_6_1.eContents().get(0);
		private final Group cGroup_8_6_2 = (Group)cGroup_8_6.eContents().get(2);
		private final Keyword cCommaKeyword_8_6_2_0 = (Keyword)cGroup_8_6_2.eContents().get(0);
		private final Assignment cSupportedFeaturesAssignment_8_6_2_1 = (Assignment)cGroup_8_6_2.eContents().get(1);
		private final RuleCall cSupportedFeaturesMdGraphFeatureEnumRuleCall_8_6_2_1_0 = (RuleCall)cSupportedFeaturesAssignment_8_6_2_1.eContents().get(0);
		private final Group cGroup_8_7 = (Group)cUnorderedGroup_8.eContents().get(7);
		private final Keyword cValidatorKeyword_8_7_0 = (Keyword)cGroup_8_7.eContents().get(0);
		private final Assignment cValidatorAssignment_8_7_1 = (Assignment)cGroup_8_7.eContents().get(1);
		private final RuleCall cValidatorJvmTypeReferenceParserRuleCall_8_7_1_0 = (RuleCall)cValidatorAssignment_8_7_1.eContents().get(0);
		private final Assignment cSupportedOptionsAssignment_9 = (Assignment)cGroup.eContents().get(9);
		private final RuleCall cSupportedOptionsMdOptionSupportParserRuleCall_9_0 = (RuleCall)cSupportedOptionsAssignment_9.eContents().get(0);
		private final Keyword cRightCurlyBracketKeyword_10 = (Keyword)cGroup.eContents().get(10);
		
		//MdAlgorithm:
		//	deprecated?='deprecated'?
		//	'algorithm' name=ID '(' provider=JvmTypeReference ('#' parameter=ID)? ')' '{' (('label' label=STRING)?
		//	& ('metadataClass' targetClass=QualifiedName)?
		//	& ('description' description=STRING)?
		//	& ('documentation' documentation=STRING)?
		//	& ('category' category=[MdCategory|QualifiedName])?
		//	& ('preview' previewImage=Path)?
		//	& ('features' supportedFeatures+=MdGraphFeature (',' supportedFeatures+=MdGraphFeature)*)?
		//	& ('validator' validator=JvmTypeReference)?) supportedOptions+=MdOptionSupport*
		//	'}';
		@Override public ParserRule getRule() { return rule; }
		
		//deprecated?='deprecated'? 'algorithm' name=ID '(' provider=JvmTypeReference ('#' parameter=ID)? ')' '{' (('label'
		//label=STRING)? & ('metadataClass' targetClass=QualifiedName)? & ('description' description=STRING)? & ('documentation'
		//documentation=STRING)? & ('category' category=[MdCategory|QualifiedName])? & ('preview' previewImage=Path)? &
		//('features' supportedFeatures+=MdGraphFeature (',' supportedFeatures+=MdGraphFeature)*)? & ('validator'
		//validator=JvmTypeReference)?) supportedOptions+=MdOptionSupport* '}'
		public Group getGroup() { return cGroup; }
		
		//deprecated?='deprecated'?
		public Assignment getDeprecatedAssignment_0() { return cDeprecatedAssignment_0; }
		
		//'deprecated'
		public Keyword getDeprecatedDeprecatedKeyword_0_0() { return cDeprecatedDeprecatedKeyword_0_0; }
		
		//'algorithm'
		public Keyword getAlgorithmKeyword_1() { return cAlgorithmKeyword_1; }
		
		//name=ID
		public Assignment getNameAssignment_2() { return cNameAssignment_2; }
		
		//ID
		public RuleCall getNameIDTerminalRuleCall_2_0() { return cNameIDTerminalRuleCall_2_0; }
		
		//'('
		public Keyword getLeftParenthesisKeyword_3() { return cLeftParenthesisKeyword_3; }
		
		//provider=JvmTypeReference
		public Assignment getProviderAssignment_4() { return cProviderAssignment_4; }
		
		//JvmTypeReference
		public RuleCall getProviderJvmTypeReferenceParserRuleCall_4_0() { return cProviderJvmTypeReferenceParserRuleCall_4_0; }
		
		//('#' parameter=ID)?
		public Group getGroup_5() { return cGroup_5; }
		
		//'#'
		public Keyword getNumberSignKeyword_5_0() { return cNumberSignKeyword_5_0; }
		
		//parameter=ID
		public Assignment getParameterAssignment_5_1() { return cParameterAssignment_5_1; }
		
		//ID
		public RuleCall getParameterIDTerminalRuleCall_5_1_0() { return cParameterIDTerminalRuleCall_5_1_0; }
		
		//')'
		public Keyword getRightParenthesisKeyword_6() { return cRightParenthesisKeyword_6; }
		
		//'{'
		public Keyword getLeftCurlyBracketKeyword_7() { return cLeftCurlyBracketKeyword_7; }
		
		//('label' label=STRING)? & ('metadataClass' targetClass=QualifiedName)? & ('description' description=STRING)? &
		//('documentation' documentation=STRING)? & ('category' category=[MdCategory|QualifiedName])? & ('preview'
		//previewImage=Path)? & ('features' supportedFeatures+=MdGraphFeature (',' supportedFeatures+=MdGraphFeature)*)? &
		//('validator' validator=JvmTypeReference)?
		public UnorderedGroup getUnorderedGroup_8() { return cUnorderedGroup_8; }
		
		//('label' label=STRING)?
		public Group getGroup_8_0() { return cGroup_8_0; }
		
		//'label'
		public Keyword getLabelKeyword_8_0_0() { return cLabelKeyword_8_0_0; }
		
		//label=STRING
		public Assignment getLabelAssignment_8_0_1() { return cLabelAssignment_8_0_1; }
		
		//STRING
		public RuleCall getLabelSTRINGTerminalRuleCall_8_0_1_0() { return cLabelSTRINGTerminalRuleCall_8_0_1_0; }
		
		//('metadataClass' targetClass=QualifiedName)?
		public Group getGroup_8_1() { return cGroup_8_1; }
		
		//'metadataClass'
		public Keyword getMetadataClassKeyword_8_1_0() { return cMetadataClassKeyword_8_1_0; }
		
		//targetClass=QualifiedName
		public Assignment getTargetClassAssignment_8_1_1() { return cTargetClassAssignment_8_1_1; }
		
		//QualifiedName
		public RuleCall getTargetClassQualifiedNameParserRuleCall_8_1_1_0() { return cTargetClassQualifiedNameParserRuleCall_8_1_1_0; }
		
		//('description' description=STRING)?
		public Group getGroup_8_2() { return cGroup_8_2; }
		
		//'description'
		public Keyword getDescriptionKeyword_8_2_0() { return cDescriptionKeyword_8_2_0; }
		
		//description=STRING
		public Assignment getDescriptionAssignment_8_2_1() { return cDescriptionAssignment_8_2_1; }
		
		//STRING
		public RuleCall getDescriptionSTRINGTerminalRuleCall_8_2_1_0() { return cDescriptionSTRINGTerminalRuleCall_8_2_1_0; }
		
		//('documentation' documentation=STRING)?
		public Group getGroup_8_3() { return cGroup_8_3; }
		
		//'documentation'
		public Keyword getDocumentationKeyword_8_3_0() { return cDocumentationKeyword_8_3_0; }
		
		//documentation=STRING
		public Assignment getDocumentationAssignment_8_3_1() { return cDocumentationAssignment_8_3_1; }
		
		//STRING
		public RuleCall getDocumentationSTRINGTerminalRuleCall_8_3_1_0() { return cDocumentationSTRINGTerminalRuleCall_8_3_1_0; }
		
		//('category' category=[MdCategory|QualifiedName])?
		public Group getGroup_8_4() { return cGroup_8_4; }
		
		//'category'
		public Keyword getCategoryKeyword_8_4_0() { return cCategoryKeyword_8_4_0; }
		
		//category=[MdCategory|QualifiedName]
		public Assignment getCategoryAssignment_8_4_1() { return cCategoryAssignment_8_4_1; }
		
		//[MdCategory|QualifiedName]
		public CrossReference getCategoryMdCategoryCrossReference_8_4_1_0() { return cCategoryMdCategoryCrossReference_8_4_1_0; }
		
		//QualifiedName
		public RuleCall getCategoryMdCategoryQualifiedNameParserRuleCall_8_4_1_0_1() { return cCategoryMdCategoryQualifiedNameParserRuleCall_8_4_1_0_1; }
		
		//('preview' previewImage=Path)?
		public Group getGroup_8_5() { return cGroup_8_5; }
		
		//'preview'
		public Keyword getPreviewKeyword_8_5_0() { return cPreviewKeyword_8_5_0; }
		
		//previewImage=Path
		public Assignment getPreviewImageAssignment_8_5_1() { return cPreviewImageAssignment_8_5_1; }
		
		//Path
		public RuleCall getPreviewImagePathParserRuleCall_8_5_1_0() { return cPreviewImagePathParserRuleCall_8_5_1_0; }
		
		//('features' supportedFeatures+=MdGraphFeature (',' supportedFeatures+=MdGraphFeature)*)?
		public Group getGroup_8_6() { return cGroup_8_6; }
		
		//'features'
		public Keyword getFeaturesKeyword_8_6_0() { return cFeaturesKeyword_8_6_0; }
		
		//supportedFeatures+=MdGraphFeature
		public Assignment getSupportedFeaturesAssignment_8_6_1() { return cSupportedFeaturesAssignment_8_6_1; }
		
		//MdGraphFeature
		public RuleCall getSupportedFeaturesMdGraphFeatureEnumRuleCall_8_6_1_0() { return cSupportedFeaturesMdGraphFeatureEnumRuleCall_8_6_1_0; }
		
		//(',' supportedFeatures+=MdGraphFeature)*
		public Group getGroup_8_6_2() { return cGroup_8_6_2; }
		
		//','
		public Keyword getCommaKeyword_8_6_2_0() { return cCommaKeyword_8_6_2_0; }
		
		//supportedFeatures+=MdGraphFeature
		public Assignment getSupportedFeaturesAssignment_8_6_2_1() { return cSupportedFeaturesAssignment_8_6_2_1; }
		
		//MdGraphFeature
		public RuleCall getSupportedFeaturesMdGraphFeatureEnumRuleCall_8_6_2_1_0() { return cSupportedFeaturesMdGraphFeatureEnumRuleCall_8_6_2_1_0; }
		
		//('validator' validator=JvmTypeReference)?
		public Group getGroup_8_7() { return cGroup_8_7; }
		
		//'validator'
		public Keyword getValidatorKeyword_8_7_0() { return cValidatorKeyword_8_7_0; }
		
		//validator=JvmTypeReference
		public Assignment getValidatorAssignment_8_7_1() { return cValidatorAssignment_8_7_1; }
		
		//JvmTypeReference
		public RuleCall getValidatorJvmTypeReferenceParserRuleCall_8_7_1_0() { return cValidatorJvmTypeReferenceParserRuleCall_8_7_1_0; }
		
		//supportedOptions+=MdOptionSupport*
		public Assignment getSupportedOptionsAssignment_9() { return cSupportedOptionsAssignment_9; }
		
		//MdOptionSupport
		public RuleCall getSupportedOptionsMdOptionSupportParserRuleCall_9_0() { return cSupportedOptionsMdOptionSupportParserRuleCall_9_0; }
		
		//'}'
		public Keyword getRightCurlyBracketKeyword_10() { return cRightCurlyBracketKeyword_10; }
	}
	public class MdCategoryElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.core.meta.MetaData.MdCategory");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Assignment cDeprecatedAssignment_0 = (Assignment)cGroup.eContents().get(0);
		private final Keyword cDeprecatedDeprecatedKeyword_0_0 = (Keyword)cDeprecatedAssignment_0.eContents().get(0);
		private final Keyword cCategoryKeyword_1 = (Keyword)cGroup.eContents().get(1);
		private final Assignment cNameAssignment_2 = (Assignment)cGroup.eContents().get(2);
		private final RuleCall cNameIDTerminalRuleCall_2_0 = (RuleCall)cNameAssignment_2.eContents().get(0);
		private final Keyword cLeftCurlyBracketKeyword_3 = (Keyword)cGroup.eContents().get(3);
		private final UnorderedGroup cUnorderedGroup_4 = (UnorderedGroup)cGroup.eContents().get(4);
		private final Group cGroup_4_0 = (Group)cUnorderedGroup_4.eContents().get(0);
		private final Keyword cLabelKeyword_4_0_0 = (Keyword)cGroup_4_0.eContents().get(0);
		private final Assignment cLabelAssignment_4_0_1 = (Assignment)cGroup_4_0.eContents().get(1);
		private final RuleCall cLabelSTRINGTerminalRuleCall_4_0_1_0 = (RuleCall)cLabelAssignment_4_0_1.eContents().get(0);
		private final Group cGroup_4_1 = (Group)cUnorderedGroup_4.eContents().get(1);
		private final Keyword cDescriptionKeyword_4_1_0 = (Keyword)cGroup_4_1.eContents().get(0);
		private final Assignment cDescriptionAssignment_4_1_1 = (Assignment)cGroup_4_1.eContents().get(1);
		private final RuleCall cDescriptionSTRINGTerminalRuleCall_4_1_1_0 = (RuleCall)cDescriptionAssignment_4_1_1.eContents().get(0);
		private final Group cGroup_4_2 = (Group)cUnorderedGroup_4.eContents().get(2);
		private final Keyword cDocumentationKeyword_4_2_0 = (Keyword)cGroup_4_2.eContents().get(0);
		private final Assignment cDocumentationAssignment_4_2_1 = (Assignment)cGroup_4_2.eContents().get(1);
		private final RuleCall cDocumentationSTRINGTerminalRuleCall_4_2_1_0 = (RuleCall)cDocumentationAssignment_4_2_1.eContents().get(0);
		private final Keyword cRightCurlyBracketKeyword_5 = (Keyword)cGroup.eContents().get(5);
		
		//MdCategory:
		//	deprecated?='deprecated'?
		//	'category' name=ID '{' (('label' label=STRING)?
		//	& ('description' description=STRING)?
		//	& ('documentation' documentation=STRING)?)
		//	'}';
		@Override public ParserRule getRule() { return rule; }
		
		//deprecated?='deprecated'? 'category' name=ID '{' (('label' label=STRING)? & ('description' description=STRING)? &
		//('documentation' documentation=STRING)?) '}'
		public Group getGroup() { return cGroup; }
		
		//deprecated?='deprecated'?
		public Assignment getDeprecatedAssignment_0() { return cDeprecatedAssignment_0; }
		
		//'deprecated'
		public Keyword getDeprecatedDeprecatedKeyword_0_0() { return cDeprecatedDeprecatedKeyword_0_0; }
		
		//'category'
		public Keyword getCategoryKeyword_1() { return cCategoryKeyword_1; }
		
		//name=ID
		public Assignment getNameAssignment_2() { return cNameAssignment_2; }
		
		//ID
		public RuleCall getNameIDTerminalRuleCall_2_0() { return cNameIDTerminalRuleCall_2_0; }
		
		//'{'
		public Keyword getLeftCurlyBracketKeyword_3() { return cLeftCurlyBracketKeyword_3; }
		
		//('label' label=STRING)? & ('description' description=STRING)? & ('documentation' documentation=STRING)?
		public UnorderedGroup getUnorderedGroup_4() { return cUnorderedGroup_4; }
		
		//('label' label=STRING)?
		public Group getGroup_4_0() { return cGroup_4_0; }
		
		//'label'
		public Keyword getLabelKeyword_4_0_0() { return cLabelKeyword_4_0_0; }
		
		//label=STRING
		public Assignment getLabelAssignment_4_0_1() { return cLabelAssignment_4_0_1; }
		
		//STRING
		public RuleCall getLabelSTRINGTerminalRuleCall_4_0_1_0() { return cLabelSTRINGTerminalRuleCall_4_0_1_0; }
		
		//('description' description=STRING)?
		public Group getGroup_4_1() { return cGroup_4_1; }
		
		//'description'
		public Keyword getDescriptionKeyword_4_1_0() { return cDescriptionKeyword_4_1_0; }
		
		//description=STRING
		public Assignment getDescriptionAssignment_4_1_1() { return cDescriptionAssignment_4_1_1; }
		
		//STRING
		public RuleCall getDescriptionSTRINGTerminalRuleCall_4_1_1_0() { return cDescriptionSTRINGTerminalRuleCall_4_1_1_0; }
		
		//('documentation' documentation=STRING)?
		public Group getGroup_4_2() { return cGroup_4_2; }
		
		//'documentation'
		public Keyword getDocumentationKeyword_4_2_0() { return cDocumentationKeyword_4_2_0; }
		
		//documentation=STRING
		public Assignment getDocumentationAssignment_4_2_1() { return cDocumentationAssignment_4_2_1; }
		
		//STRING
		public RuleCall getDocumentationSTRINGTerminalRuleCall_4_2_1_0() { return cDocumentationSTRINGTerminalRuleCall_4_2_1_0; }
		
		//'}'
		public Keyword getRightCurlyBracketKeyword_5() { return cRightCurlyBracketKeyword_5; }
	}
	public class MdOptionSupportElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.core.meta.MetaData.MdOptionSupport");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final Keyword cSupportsKeyword_0 = (Keyword)cGroup.eContents().get(0);
		private final Assignment cOptionAssignment_1 = (Assignment)cGroup.eContents().get(1);
		private final CrossReference cOptionMdOptionCrossReference_1_0 = (CrossReference)cOptionAssignment_1.eContents().get(0);
		private final RuleCall cOptionMdOptionQualifiedNameParserRuleCall_1_0_1 = (RuleCall)cOptionMdOptionCrossReference_1_0.eContents().get(1);
		private final Group cGroup_2 = (Group)cGroup.eContents().get(2);
		private final Keyword cEqualsSignKeyword_2_0 = (Keyword)cGroup_2.eContents().get(0);
		private final Assignment cValueAssignment_2_1 = (Assignment)cGroup_2.eContents().get(1);
		private final RuleCall cValueXExpressionParserRuleCall_2_1_0 = (RuleCall)cValueAssignment_2_1.eContents().get(0);
		private final Group cGroup_3 = (Group)cGroup.eContents().get(3);
		private final Keyword cDocumentationKeyword_3_0 = (Keyword)cGroup_3.eContents().get(0);
		private final Assignment cDocumentationAssignment_3_1 = (Assignment)cGroup_3.eContents().get(1);
		private final RuleCall cDocumentationSTRINGTerminalRuleCall_3_1_0 = (RuleCall)cDocumentationAssignment_3_1.eContents().get(0);
		
		//MdOptionSupport:
		//	'supports' option=[MdOption|QualifiedName] ('=' value=XExpression)? ('documentation' documentation=STRING)?;
		@Override public ParserRule getRule() { return rule; }
		
		//'supports' option=[MdOption|QualifiedName] ('=' value=XExpression)? ('documentation' documentation=STRING)?
		public Group getGroup() { return cGroup; }
		
		//'supports'
		public Keyword getSupportsKeyword_0() { return cSupportsKeyword_0; }
		
		//option=[MdOption|QualifiedName]
		public Assignment getOptionAssignment_1() { return cOptionAssignment_1; }
		
		//[MdOption|QualifiedName]
		public CrossReference getOptionMdOptionCrossReference_1_0() { return cOptionMdOptionCrossReference_1_0; }
		
		//QualifiedName
		public RuleCall getOptionMdOptionQualifiedNameParserRuleCall_1_0_1() { return cOptionMdOptionQualifiedNameParserRuleCall_1_0_1; }
		
		//('=' value=XExpression)?
		public Group getGroup_2() { return cGroup_2; }
		
		//'='
		public Keyword getEqualsSignKeyword_2_0() { return cEqualsSignKeyword_2_0; }
		
		//value=XExpression
		public Assignment getValueAssignment_2_1() { return cValueAssignment_2_1; }
		
		//XExpression
		public RuleCall getValueXExpressionParserRuleCall_2_1_0() { return cValueXExpressionParserRuleCall_2_1_0; }
		
		//('documentation' documentation=STRING)?
		public Group getGroup_3() { return cGroup_3; }
		
		//'documentation'
		public Keyword getDocumentationKeyword_3_0() { return cDocumentationKeyword_3_0; }
		
		//documentation=STRING
		public Assignment getDocumentationAssignment_3_1() { return cDocumentationAssignment_3_1; }
		
		//STRING
		public RuleCall getDocumentationSTRINGTerminalRuleCall_3_1_0() { return cDocumentationSTRINGTerminalRuleCall_3_1_0; }
	}
	public class PathElements extends AbstractParserRuleElementFinder {
		private final ParserRule rule = (ParserRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.core.meta.MetaData.Path");
		private final Group cGroup = (Group)rule.eContents().get(1);
		private final RuleCall cQualifiedNameParserRuleCall_0 = (RuleCall)cGroup.eContents().get(0);
		private final Group cGroup_1 = (Group)cGroup.eContents().get(1);
		private final Alternatives cAlternatives_1_0 = (Alternatives)cGroup_1.eContents().get(0);
		private final Keyword cSolidusKeyword_1_0_0 = (Keyword)cAlternatives_1_0.eContents().get(0);
		private final Keyword cHyphenMinusKeyword_1_0_1 = (Keyword)cAlternatives_1_0.eContents().get(1);
		private final RuleCall cQualifiedNameParserRuleCall_1_1 = (RuleCall)cGroup_1.eContents().get(1);
		
		//Path:
		//	QualifiedName (('/' | '-') QualifiedName)*;
		@Override public ParserRule getRule() { return rule; }
		
		//QualifiedName (('/' | '-') QualifiedName)*
		public Group getGroup() { return cGroup; }
		
		//QualifiedName
		public RuleCall getQualifiedNameParserRuleCall_0() { return cQualifiedNameParserRuleCall_0; }
		
		//(('/' | '-') QualifiedName)*
		public Group getGroup_1() { return cGroup_1; }
		
		//'/' | '-'
		public Alternatives getAlternatives_1_0() { return cAlternatives_1_0; }
		
		//'/'
		public Keyword getSolidusKeyword_1_0_0() { return cSolidusKeyword_1_0_0; }
		
		//'-'
		public Keyword getHyphenMinusKeyword_1_0_1() { return cHyphenMinusKeyword_1_0_1; }
		
		//QualifiedName
		public RuleCall getQualifiedNameParserRuleCall_1_1() { return cQualifiedNameParserRuleCall_1_1; }
	}
	
	public class MdOptionTargetTypeElements extends AbstractEnumRuleElementFinder {
		private final EnumRule rule = (EnumRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.core.meta.MetaData.MdOptionTargetType");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final EnumLiteralDeclaration cParentsEnumLiteralDeclaration_0 = (EnumLiteralDeclaration)cAlternatives.eContents().get(0);
		private final Keyword cParentsParentsKeyword_0_0 = (Keyword)cParentsEnumLiteralDeclaration_0.eContents().get(0);
		private final EnumLiteralDeclaration cNodesEnumLiteralDeclaration_1 = (EnumLiteralDeclaration)cAlternatives.eContents().get(1);
		private final Keyword cNodesNodesKeyword_1_0 = (Keyword)cNodesEnumLiteralDeclaration_1.eContents().get(0);
		private final EnumLiteralDeclaration cEdgesEnumLiteralDeclaration_2 = (EnumLiteralDeclaration)cAlternatives.eContents().get(2);
		private final Keyword cEdgesEdgesKeyword_2_0 = (Keyword)cEdgesEnumLiteralDeclaration_2.eContents().get(0);
		private final EnumLiteralDeclaration cPortsEnumLiteralDeclaration_3 = (EnumLiteralDeclaration)cAlternatives.eContents().get(3);
		private final Keyword cPortsPortsKeyword_3_0 = (Keyword)cPortsEnumLiteralDeclaration_3.eContents().get(0);
		private final EnumLiteralDeclaration cLabelsEnumLiteralDeclaration_4 = (EnumLiteralDeclaration)cAlternatives.eContents().get(4);
		private final Keyword cLabelsLabelsKeyword_4_0 = (Keyword)cLabelsEnumLiteralDeclaration_4.eContents().get(0);
		
		//enum MdOptionTargetType:
		//	parents | nodes | edges | ports | labels;
		public EnumRule getRule() { return rule; }
		
		//parents | nodes | edges | ports | labels
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//parents
		public EnumLiteralDeclaration getParentsEnumLiteralDeclaration_0() { return cParentsEnumLiteralDeclaration_0; }
		
		//"parents"
		public Keyword getParentsParentsKeyword_0_0() { return cParentsParentsKeyword_0_0; }
		
		//nodes
		public EnumLiteralDeclaration getNodesEnumLiteralDeclaration_1() { return cNodesEnumLiteralDeclaration_1; }
		
		//"nodes"
		public Keyword getNodesNodesKeyword_1_0() { return cNodesNodesKeyword_1_0; }
		
		//edges
		public EnumLiteralDeclaration getEdgesEnumLiteralDeclaration_2() { return cEdgesEnumLiteralDeclaration_2; }
		
		//"edges"
		public Keyword getEdgesEdgesKeyword_2_0() { return cEdgesEdgesKeyword_2_0; }
		
		//ports
		public EnumLiteralDeclaration getPortsEnumLiteralDeclaration_3() { return cPortsEnumLiteralDeclaration_3; }
		
		//"ports"
		public Keyword getPortsPortsKeyword_3_0() { return cPortsPortsKeyword_3_0; }
		
		//labels
		public EnumLiteralDeclaration getLabelsEnumLiteralDeclaration_4() { return cLabelsEnumLiteralDeclaration_4; }
		
		//"labels"
		public Keyword getLabelsLabelsKeyword_4_0() { return cLabelsLabelsKeyword_4_0; }
	}
	public class MdGraphFeatureElements extends AbstractEnumRuleElementFinder {
		private final EnumRule rule = (EnumRule) GrammarUtil.findRuleForName(getGrammar(), "org.eclipse.elk.core.meta.MetaData.MdGraphFeature");
		private final Alternatives cAlternatives = (Alternatives)rule.eContents().get(1);
		private final EnumLiteralDeclaration cSelf_loopsEnumLiteralDeclaration_0 = (EnumLiteralDeclaration)cAlternatives.eContents().get(0);
		private final Keyword cSelf_loopsSelf_loopsKeyword_0_0 = (Keyword)cSelf_loopsEnumLiteralDeclaration_0.eContents().get(0);
		private final EnumLiteralDeclaration cInside_self_loopsEnumLiteralDeclaration_1 = (EnumLiteralDeclaration)cAlternatives.eContents().get(1);
		private final Keyword cInside_self_loopsInside_self_loopsKeyword_1_0 = (Keyword)cInside_self_loopsEnumLiteralDeclaration_1.eContents().get(0);
		private final EnumLiteralDeclaration cMulti_edgesEnumLiteralDeclaration_2 = (EnumLiteralDeclaration)cAlternatives.eContents().get(2);
		private final Keyword cMulti_edgesMulti_edgesKeyword_2_0 = (Keyword)cMulti_edgesEnumLiteralDeclaration_2.eContents().get(0);
		private final EnumLiteralDeclaration cEdge_labelsEnumLiteralDeclaration_3 = (EnumLiteralDeclaration)cAlternatives.eContents().get(3);
		private final Keyword cEdge_labelsEdge_labelsKeyword_3_0 = (Keyword)cEdge_labelsEnumLiteralDeclaration_3.eContents().get(0);
		private final EnumLiteralDeclaration cPortsEnumLiteralDeclaration_4 = (EnumLiteralDeclaration)cAlternatives.eContents().get(4);
		private final Keyword cPortsPortsKeyword_4_0 = (Keyword)cPortsEnumLiteralDeclaration_4.eContents().get(0);
		private final EnumLiteralDeclaration cCompoundEnumLiteralDeclaration_5 = (EnumLiteralDeclaration)cAlternatives.eContents().get(5);
		private final Keyword cCompoundCompoundKeyword_5_0 = (Keyword)cCompoundEnumLiteralDeclaration_5.eContents().get(0);
		private final EnumLiteralDeclaration cClustersEnumLiteralDeclaration_6 = (EnumLiteralDeclaration)cAlternatives.eContents().get(6);
		private final Keyword cClustersClustersKeyword_6_0 = (Keyword)cClustersEnumLiteralDeclaration_6.eContents().get(0);
		private final EnumLiteralDeclaration cDisconnectedEnumLiteralDeclaration_7 = (EnumLiteralDeclaration)cAlternatives.eContents().get(7);
		private final Keyword cDisconnectedDisconnectedKeyword_7_0 = (Keyword)cDisconnectedEnumLiteralDeclaration_7.eContents().get(0);
		
		//enum MdGraphFeature:
		//	self_loops | inside_self_loops | multi_edges | edge_labels | ports | compound
		//	| clusters | disconnected;
		public EnumRule getRule() { return rule; }
		
		//self_loops | inside_self_loops | multi_edges | edge_labels | ports | compound | clusters | disconnected
		public Alternatives getAlternatives() { return cAlternatives; }
		
		//self_loops
		public EnumLiteralDeclaration getSelf_loopsEnumLiteralDeclaration_0() { return cSelf_loopsEnumLiteralDeclaration_0; }
		
		//"self_loops"
		public Keyword getSelf_loopsSelf_loopsKeyword_0_0() { return cSelf_loopsSelf_loopsKeyword_0_0; }
		
		//inside_self_loops
		public EnumLiteralDeclaration getInside_self_loopsEnumLiteralDeclaration_1() { return cInside_self_loopsEnumLiteralDeclaration_1; }
		
		//"inside_self_loops"
		public Keyword getInside_self_loopsInside_self_loopsKeyword_1_0() { return cInside_self_loopsInside_self_loopsKeyword_1_0; }
		
		//multi_edges
		public EnumLiteralDeclaration getMulti_edgesEnumLiteralDeclaration_2() { return cMulti_edgesEnumLiteralDeclaration_2; }
		
		//"multi_edges"
		public Keyword getMulti_edgesMulti_edgesKeyword_2_0() { return cMulti_edgesMulti_edgesKeyword_2_0; }
		
		//edge_labels
		public EnumLiteralDeclaration getEdge_labelsEnumLiteralDeclaration_3() { return cEdge_labelsEnumLiteralDeclaration_3; }
		
		//"edge_labels"
		public Keyword getEdge_labelsEdge_labelsKeyword_3_0() { return cEdge_labelsEdge_labelsKeyword_3_0; }
		
		//ports
		public EnumLiteralDeclaration getPortsEnumLiteralDeclaration_4() { return cPortsEnumLiteralDeclaration_4; }
		
		//"ports"
		public Keyword getPortsPortsKeyword_4_0() { return cPortsPortsKeyword_4_0; }
		
		//compound
		public EnumLiteralDeclaration getCompoundEnumLiteralDeclaration_5() { return cCompoundEnumLiteralDeclaration_5; }
		
		//"compound"
		public Keyword getCompoundCompoundKeyword_5_0() { return cCompoundCompoundKeyword_5_0; }
		
		//clusters
		public EnumLiteralDeclaration getClustersEnumLiteralDeclaration_6() { return cClustersEnumLiteralDeclaration_6; }
		
		//"clusters"
		public Keyword getClustersClustersKeyword_6_0() { return cClustersClustersKeyword_6_0; }
		
		//disconnected
		public EnumLiteralDeclaration getDisconnectedEnumLiteralDeclaration_7() { return cDisconnectedEnumLiteralDeclaration_7; }
		
		//"disconnected"
		public Keyword getDisconnectedDisconnectedKeyword_7_0() { return cDisconnectedDisconnectedKeyword_7_0; }
	}
	
	private final MdModelElements pMdModel;
	private final MdBundleElements pMdBundle;
	private final MdBundleMemberElements pMdBundleMember;
	private final MdGroupOrOptionElements pMdGroupOrOption;
	private final MdGroupElements pMdGroup;
	private final MdOptionElements pMdOption;
	private final MdOptionDependencyElements pMdOptionDependency;
	private final MdAlgorithmElements pMdAlgorithm;
	private final MdCategoryElements pMdCategory;
	private final MdOptionSupportElements pMdOptionSupport;
	private final PathElements pPath;
	private final MdOptionTargetTypeElements eMdOptionTargetType;
	private final MdGraphFeatureElements eMdGraphFeature;
	
	private final Grammar grammar;
	
	private final XbaseGrammarAccess gaXbase;
	
	private final XtypeGrammarAccess gaXtype;

	@Inject
	public MetaDataGrammarAccess(GrammarProvider grammarProvider,
			XbaseGrammarAccess gaXbase,
			XtypeGrammarAccess gaXtype) {
		this.grammar = internalFindGrammar(grammarProvider);
		this.gaXbase = gaXbase;
		this.gaXtype = gaXtype;
		this.pMdModel = new MdModelElements();
		this.pMdBundle = new MdBundleElements();
		this.pMdBundleMember = new MdBundleMemberElements();
		this.pMdGroupOrOption = new MdGroupOrOptionElements();
		this.pMdGroup = new MdGroupElements();
		this.pMdOption = new MdOptionElements();
		this.pMdOptionDependency = new MdOptionDependencyElements();
		this.pMdAlgorithm = new MdAlgorithmElements();
		this.pMdCategory = new MdCategoryElements();
		this.pMdOptionSupport = new MdOptionSupportElements();
		this.pPath = new PathElements();
		this.eMdOptionTargetType = new MdOptionTargetTypeElements();
		this.eMdGraphFeature = new MdGraphFeatureElements();
	}
	
	protected Grammar internalFindGrammar(GrammarProvider grammarProvider) {
		Grammar grammar = grammarProvider.getGrammar(this);
		while (grammar != null) {
			if ("org.eclipse.elk.core.meta.MetaData".equals(grammar.getName())) {
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
	
	
	public XbaseGrammarAccess getXbaseGrammarAccess() {
		return gaXbase;
	}
	
	public XtypeGrammarAccess getXtypeGrammarAccess() {
		return gaXtype;
	}

	
	//MdModel:
	//	('package' name=QualifiedName
	//	importSection=XImportSection?
	//	bundle=MdBundle)?;
	public MdModelElements getMdModelAccess() {
		return pMdModel;
	}
	
	public ParserRule getMdModelRule() {
		return getMdModelAccess().getRule();
	}
	
	//MdBundle:
	//	{MdBundle} ('bundle' '{' (('label' label=STRING)?
	//	& ('metadataClass' targetClass=QualifiedName)?
	//	& ('documentationFolder' documentationFolder=Path)?
	//	& ('idPrefix' idPrefix=QualifiedName)?)
	//	'}')?
	//	members+=MdBundleMember*;
	public MdBundleElements getMdBundleAccess() {
		return pMdBundle;
	}
	
	public ParserRule getMdBundleRule() {
		return getMdBundleAccess().getRule();
	}
	
	//MdBundleMember:
	//	MdGroupOrOption | MdAlgorithm | MdCategory;
	public MdBundleMemberElements getMdBundleMemberAccess() {
		return pMdBundleMember;
	}
	
	public ParserRule getMdBundleMemberRule() {
		return getMdBundleMemberAccess().getRule();
	}
	
	//MdGroupOrOption:
	//	MdGroup | MdOption;
	public MdGroupOrOptionElements getMdGroupOrOptionAccess() {
		return pMdGroupOrOption;
	}
	
	public ParserRule getMdGroupOrOptionRule() {
		return getMdGroupOrOptionAccess().getRule();
	}
	
	//MdGroup:
	//	'group' name=ID '{' ('documentation' documentation=STRING)?
	//	children+=MdGroupOrOption*
	//	'}';
	public MdGroupElements getMdGroupAccess() {
		return pMdGroup;
	}
	
	public ParserRule getMdGroupRule() {
		return getMdGroupAccess().getRule();
	}
	
	//MdOption:
	//	deprecated?='deprecated'? (advanced?='advanced' | programmatic?='programmatic' | output?='output' | global?='global')?
	//	'option' name=ID (':' type=JvmTypeReference)? '{' (('label' label=STRING)?
	//	& ('description' description=STRING)?
	//	& ('documentation' documentation=STRING)? // "@docu/priority.md"
	//	& ('default' '=' defaultValue=XExpression)?
	//	& ('lowerBound' '=' lowerBound=XExpression)?
	//	& ('upperBound' '=' upperBound=XExpression)?
	//	& ('targets' targets+=MdOptionTargetType (',' targets+=MdOptionTargetType)*)?
	//	& ('legacyIds' legacyIds+=QualifiedName (',' legacyIds+=QualifiedName)*)?) dependencies+=MdOptionDependency*
	//	'}';
	public MdOptionElements getMdOptionAccess() {
		return pMdOption;
	}
	
	public ParserRule getMdOptionRule() {
		return getMdOptionAccess().getRule();
	}
	
	//MdOptionDependency:
	//	'requires' target=[MdOption|QualifiedName] ('==' value=XExpression)?;
	public MdOptionDependencyElements getMdOptionDependencyAccess() {
		return pMdOptionDependency;
	}
	
	public ParserRule getMdOptionDependencyRule() {
		return getMdOptionDependencyAccess().getRule();
	}
	
	//MdAlgorithm:
	//	deprecated?='deprecated'?
	//	'algorithm' name=ID '(' provider=JvmTypeReference ('#' parameter=ID)? ')' '{' (('label' label=STRING)?
	//	& ('metadataClass' targetClass=QualifiedName)?
	//	& ('description' description=STRING)?
	//	& ('documentation' documentation=STRING)?
	//	& ('category' category=[MdCategory|QualifiedName])?
	//	& ('preview' previewImage=Path)?
	//	& ('features' supportedFeatures+=MdGraphFeature (',' supportedFeatures+=MdGraphFeature)*)?
	//	& ('validator' validator=JvmTypeReference)?) supportedOptions+=MdOptionSupport*
	//	'}';
	public MdAlgorithmElements getMdAlgorithmAccess() {
		return pMdAlgorithm;
	}
	
	public ParserRule getMdAlgorithmRule() {
		return getMdAlgorithmAccess().getRule();
	}
	
	//MdCategory:
	//	deprecated?='deprecated'?
	//	'category' name=ID '{' (('label' label=STRING)?
	//	& ('description' description=STRING)?
	//	& ('documentation' documentation=STRING)?)
	//	'}';
	public MdCategoryElements getMdCategoryAccess() {
		return pMdCategory;
	}
	
	public ParserRule getMdCategoryRule() {
		return getMdCategoryAccess().getRule();
	}
	
	//MdOptionSupport:
	//	'supports' option=[MdOption|QualifiedName] ('=' value=XExpression)? ('documentation' documentation=STRING)?;
	public MdOptionSupportElements getMdOptionSupportAccess() {
		return pMdOptionSupport;
	}
	
	public ParserRule getMdOptionSupportRule() {
		return getMdOptionSupportAccess().getRule();
	}
	
	//Path:
	//	QualifiedName (('/' | '-') QualifiedName)*;
	public PathElements getPathAccess() {
		return pPath;
	}
	
	public ParserRule getPathRule() {
		return getPathAccess().getRule();
	}
	
	//enum MdOptionTargetType:
	//	parents | nodes | edges | ports | labels;
	public MdOptionTargetTypeElements getMdOptionTargetTypeAccess() {
		return eMdOptionTargetType;
	}
	
	public EnumRule getMdOptionTargetTypeRule() {
		return getMdOptionTargetTypeAccess().getRule();
	}
	
	//enum MdGraphFeature:
	//	self_loops | inside_self_loops | multi_edges | edge_labels | ports | compound
	//	| clusters | disconnected;
	public MdGraphFeatureElements getMdGraphFeatureAccess() {
		return eMdGraphFeature;
	}
	
	public EnumRule getMdGraphFeatureRule() {
		return getMdGraphFeatureAccess().getRule();
	}
	
	//XExpression:
	//	XAssignment;
	public XbaseGrammarAccess.XExpressionElements getXExpressionAccess() {
		return gaXbase.getXExpressionAccess();
	}
	
	public ParserRule getXExpressionRule() {
		return getXExpressionAccess().getRule();
	}
	
	//XAssignment XExpression:
	//	{XAssignment} feature=[types::JvmIdentifiableElement|FeatureCallID] OpSingleAssign value=XAssignment | XOrExpression
	//	(=> ({XBinaryOperation.leftOperand=current} feature=[types::JvmIdentifiableElement|OpMultiAssign])
	//	rightOperand=XAssignment)?;
	public XbaseGrammarAccess.XAssignmentElements getXAssignmentAccess() {
		return gaXbase.getXAssignmentAccess();
	}
	
	public ParserRule getXAssignmentRule() {
		return getXAssignmentAccess().getRule();
	}
	
	//OpSingleAssign:
	//	'=';
	public XbaseGrammarAccess.OpSingleAssignElements getOpSingleAssignAccess() {
		return gaXbase.getOpSingleAssignAccess();
	}
	
	public ParserRule getOpSingleAssignRule() {
		return getOpSingleAssignAccess().getRule();
	}
	
	//OpMultiAssign:
	//	'+=' | '-=' | '*=' | '/=' | '%=' |
	//	'<' '<' '=' |
	//	'>' '>'? '>=';
	public XbaseGrammarAccess.OpMultiAssignElements getOpMultiAssignAccess() {
		return gaXbase.getOpMultiAssignAccess();
	}
	
	public ParserRule getOpMultiAssignRule() {
		return getOpMultiAssignAccess().getRule();
	}
	
	//XOrExpression XExpression:
	//	XAndExpression (=> ({XBinaryOperation.leftOperand=current} feature=[types::JvmIdentifiableElement|OpOr])
	//	rightOperand=XAndExpression)*;
	public XbaseGrammarAccess.XOrExpressionElements getXOrExpressionAccess() {
		return gaXbase.getXOrExpressionAccess();
	}
	
	public ParserRule getXOrExpressionRule() {
		return getXOrExpressionAccess().getRule();
	}
	
	//OpOr:
	//	'||';
	public XbaseGrammarAccess.OpOrElements getOpOrAccess() {
		return gaXbase.getOpOrAccess();
	}
	
	public ParserRule getOpOrRule() {
		return getOpOrAccess().getRule();
	}
	
	//XAndExpression XExpression:
	//	XEqualityExpression (=> ({XBinaryOperation.leftOperand=current} feature=[types::JvmIdentifiableElement|OpAnd])
	//	rightOperand=XEqualityExpression)*;
	public XbaseGrammarAccess.XAndExpressionElements getXAndExpressionAccess() {
		return gaXbase.getXAndExpressionAccess();
	}
	
	public ParserRule getXAndExpressionRule() {
		return getXAndExpressionAccess().getRule();
	}
	
	//OpAnd:
	//	'&&';
	public XbaseGrammarAccess.OpAndElements getOpAndAccess() {
		return gaXbase.getOpAndAccess();
	}
	
	public ParserRule getOpAndRule() {
		return getOpAndAccess().getRule();
	}
	
	//XEqualityExpression XExpression:
	//	XRelationalExpression (=> ({XBinaryOperation.leftOperand=current} feature=[types::JvmIdentifiableElement|OpEquality])
	//	rightOperand=XRelationalExpression)*;
	public XbaseGrammarAccess.XEqualityExpressionElements getXEqualityExpressionAccess() {
		return gaXbase.getXEqualityExpressionAccess();
	}
	
	public ParserRule getXEqualityExpressionRule() {
		return getXEqualityExpressionAccess().getRule();
	}
	
	//OpEquality:
	//	'==' | '!=' | '===' | '!==';
	public XbaseGrammarAccess.OpEqualityElements getOpEqualityAccess() {
		return gaXbase.getOpEqualityAccess();
	}
	
	public ParserRule getOpEqualityRule() {
		return getOpEqualityAccess().getRule();
	}
	
	//XRelationalExpression XExpression:
	//	XOtherOperatorExpression (=> ({XInstanceOfExpression.expression=current} 'instanceof') type=JvmTypeReference |
	//	=> ({XBinaryOperation.leftOperand=current} feature=[types::JvmIdentifiableElement|OpCompare])
	//	rightOperand=XOtherOperatorExpression)*;
	public XbaseGrammarAccess.XRelationalExpressionElements getXRelationalExpressionAccess() {
		return gaXbase.getXRelationalExpressionAccess();
	}
	
	public ParserRule getXRelationalExpressionRule() {
		return getXRelationalExpressionAccess().getRule();
	}
	
	//OpCompare:
	//	'>=' | '<' '=' | '>' | '<';
	public XbaseGrammarAccess.OpCompareElements getOpCompareAccess() {
		return gaXbase.getOpCompareAccess();
	}
	
	public ParserRule getOpCompareRule() {
		return getOpCompareAccess().getRule();
	}
	
	//XOtherOperatorExpression XExpression:
	//	XAdditiveExpression (=> ({XBinaryOperation.leftOperand=current} feature=[types::JvmIdentifiableElement|OpOther])
	//	rightOperand=XAdditiveExpression)*;
	public XbaseGrammarAccess.XOtherOperatorExpressionElements getXOtherOperatorExpressionAccess() {
		return gaXbase.getXOtherOperatorExpressionAccess();
	}
	
	public ParserRule getXOtherOperatorExpressionRule() {
		return getXOtherOperatorExpressionAccess().getRule();
	}
	
	//OpOther:
	//	'->'
	//	| '..<'
	//	| '>' '..'
	//	| '..'
	//	| '=>'
	//	| '>' (=> ('>' '>') | '>') | '<' (=> ('<' '<') | '<' | '=>') | '<>'
	//	| '?:';
	public XbaseGrammarAccess.OpOtherElements getOpOtherAccess() {
		return gaXbase.getOpOtherAccess();
	}
	
	public ParserRule getOpOtherRule() {
		return getOpOtherAccess().getRule();
	}
	
	//XAdditiveExpression XExpression:
	//	XMultiplicativeExpression (=> ({XBinaryOperation.leftOperand=current} feature=[types::JvmIdentifiableElement|OpAdd])
	//	rightOperand=XMultiplicativeExpression)*;
	public XbaseGrammarAccess.XAdditiveExpressionElements getXAdditiveExpressionAccess() {
		return gaXbase.getXAdditiveExpressionAccess();
	}
	
	public ParserRule getXAdditiveExpressionRule() {
		return getXAdditiveExpressionAccess().getRule();
	}
	
	//OpAdd:
	//	'+' | '-';
	public XbaseGrammarAccess.OpAddElements getOpAddAccess() {
		return gaXbase.getOpAddAccess();
	}
	
	public ParserRule getOpAddRule() {
		return getOpAddAccess().getRule();
	}
	
	//XMultiplicativeExpression XExpression:
	//	XUnaryOperation (=> ({XBinaryOperation.leftOperand=current} feature=[types::JvmIdentifiableElement|OpMulti])
	//	rightOperand=XUnaryOperation)*;
	public XbaseGrammarAccess.XMultiplicativeExpressionElements getXMultiplicativeExpressionAccess() {
		return gaXbase.getXMultiplicativeExpressionAccess();
	}
	
	public ParserRule getXMultiplicativeExpressionRule() {
		return getXMultiplicativeExpressionAccess().getRule();
	}
	
	//OpMulti:
	//	'*' | '**' | '/' | '%';
	public XbaseGrammarAccess.OpMultiElements getOpMultiAccess() {
		return gaXbase.getOpMultiAccess();
	}
	
	public ParserRule getOpMultiRule() {
		return getOpMultiAccess().getRule();
	}
	
	//XUnaryOperation XExpression:
	//	{XUnaryOperation} feature=[types::JvmIdentifiableElement|OpUnary] operand=XUnaryOperation
	//	| XCastedExpression;
	public XbaseGrammarAccess.XUnaryOperationElements getXUnaryOperationAccess() {
		return gaXbase.getXUnaryOperationAccess();
	}
	
	public ParserRule getXUnaryOperationRule() {
		return getXUnaryOperationAccess().getRule();
	}
	
	//OpUnary:
	//	"!" | "-" | "+";
	public XbaseGrammarAccess.OpUnaryElements getOpUnaryAccess() {
		return gaXbase.getOpUnaryAccess();
	}
	
	public ParserRule getOpUnaryRule() {
		return getOpUnaryAccess().getRule();
	}
	
	//XCastedExpression XExpression:
	//	XPostfixOperation (=> ({XCastedExpression.target=current} 'as') type=JvmTypeReference)*;
	public XbaseGrammarAccess.XCastedExpressionElements getXCastedExpressionAccess() {
		return gaXbase.getXCastedExpressionAccess();
	}
	
	public ParserRule getXCastedExpressionRule() {
		return getXCastedExpressionAccess().getRule();
	}
	
	//XPostfixOperation XExpression:
	//	XMemberFeatureCall => ({XPostfixOperation.operand=current} feature=[types::JvmIdentifiableElement|OpPostfix])?;
	public XbaseGrammarAccess.XPostfixOperationElements getXPostfixOperationAccess() {
		return gaXbase.getXPostfixOperationAccess();
	}
	
	public ParserRule getXPostfixOperationRule() {
		return getXPostfixOperationAccess().getRule();
	}
	
	//OpPostfix:
	//	"++" | "--";
	public XbaseGrammarAccess.OpPostfixElements getOpPostfixAccess() {
		return gaXbase.getOpPostfixAccess();
	}
	
	public ParserRule getOpPostfixRule() {
		return getOpPostfixAccess().getRule();
	}
	
	//XMemberFeatureCall XExpression:
	//	XPrimaryExpression (=> ({XAssignment.assignable=current} ('.' | explicitStatic?="::")
	//	feature=[types::JvmIdentifiableElement|FeatureCallID] OpSingleAssign) value=XAssignment
	//	| => ({XMemberFeatureCall.memberCallTarget=current} ("." | nullSafe?="?." | explicitStatic?="::")) ('<'
	//	typeArguments+=JvmArgumentTypeReference (',' typeArguments+=JvmArgumentTypeReference)* '>')?
	//	feature=[types::JvmIdentifiableElement|IdOrSuper] (=> explicitOperationCall?='(' (memberCallArguments+=XShortClosure
	//	| memberCallArguments+=XExpression (',' memberCallArguments+=XExpression)*)?
	//	')')?
	//	memberCallArguments+=XClosure?)*;
	public XbaseGrammarAccess.XMemberFeatureCallElements getXMemberFeatureCallAccess() {
		return gaXbase.getXMemberFeatureCallAccess();
	}
	
	public ParserRule getXMemberFeatureCallRule() {
		return getXMemberFeatureCallAccess().getRule();
	}
	
	//XPrimaryExpression XExpression:
	//	XConstructorCall | XBlockExpression | XSwitchExpression | XSynchronizedExpression | XFeatureCall | XLiteral |
	//	XIfExpression | XForLoopExpression | XBasicForLoopExpression | XWhileExpression | XDoWhileExpression |
	//	XThrowExpression | XReturnExpression | XTryCatchFinallyExpression | XParenthesizedExpression;
	public XbaseGrammarAccess.XPrimaryExpressionElements getXPrimaryExpressionAccess() {
		return gaXbase.getXPrimaryExpressionAccess();
	}
	
	public ParserRule getXPrimaryExpressionRule() {
		return getXPrimaryExpressionAccess().getRule();
	}
	
	//XLiteral XExpression:
	//	XCollectionLiteral | XClosure | XBooleanLiteral | XNumberLiteral | XNullLiteral | XStringLiteral | XTypeLiteral;
	public XbaseGrammarAccess.XLiteralElements getXLiteralAccess() {
		return gaXbase.getXLiteralAccess();
	}
	
	public ParserRule getXLiteralRule() {
		return getXLiteralAccess().getRule();
	}
	
	//XCollectionLiteral:
	//	XSetLiteral | XListLiteral;
	public XbaseGrammarAccess.XCollectionLiteralElements getXCollectionLiteralAccess() {
		return gaXbase.getXCollectionLiteralAccess();
	}
	
	public ParserRule getXCollectionLiteralRule() {
		return getXCollectionLiteralAccess().getRule();
	}
	
	//XSetLiteral:
	//	{XSetLiteral} '#' '{' (elements+=XExpression (',' elements+=XExpression)*)? '}';
	public XbaseGrammarAccess.XSetLiteralElements getXSetLiteralAccess() {
		return gaXbase.getXSetLiteralAccess();
	}
	
	public ParserRule getXSetLiteralRule() {
		return getXSetLiteralAccess().getRule();
	}
	
	//XListLiteral:
	//	{XListLiteral} '#' '[' (elements+=XExpression (',' elements+=XExpression)*)? ']';
	public XbaseGrammarAccess.XListLiteralElements getXListLiteralAccess() {
		return gaXbase.getXListLiteralAccess();
	}
	
	public ParserRule getXListLiteralRule() {
		return getXListLiteralAccess().getRule();
	}
	
	//XClosure XExpression:
	//	=> ({XClosure}
	//	'[')
	//	=> ((declaredFormalParameters+=JvmFormalParameter (',' declaredFormalParameters+=JvmFormalParameter)*)?
	//	explicitSyntax?='|')?
	//	expression=XExpressionInClosure
	//	']';
	public XbaseGrammarAccess.XClosureElements getXClosureAccess() {
		return gaXbase.getXClosureAccess();
	}
	
	public ParserRule getXClosureRule() {
		return getXClosureAccess().getRule();
	}
	
	//XExpressionInClosure XExpression:
	//	{XBlockExpression} (expressions+=XExpressionOrVarDeclaration ';'?)*;
	public XbaseGrammarAccess.XExpressionInClosureElements getXExpressionInClosureAccess() {
		return gaXbase.getXExpressionInClosureAccess();
	}
	
	public ParserRule getXExpressionInClosureRule() {
		return getXExpressionInClosureAccess().getRule();
	}
	
	//XShortClosure XExpression:
	//	=> ({XClosure} (declaredFormalParameters+=JvmFormalParameter (',' declaredFormalParameters+=JvmFormalParameter)*)?
	//	explicitSyntax?='|') expression=XExpression;
	public XbaseGrammarAccess.XShortClosureElements getXShortClosureAccess() {
		return gaXbase.getXShortClosureAccess();
	}
	
	public ParserRule getXShortClosureRule() {
		return getXShortClosureAccess().getRule();
	}
	
	//XParenthesizedExpression XExpression:
	//	'(' XExpression ')';
	public XbaseGrammarAccess.XParenthesizedExpressionElements getXParenthesizedExpressionAccess() {
		return gaXbase.getXParenthesizedExpressionAccess();
	}
	
	public ParserRule getXParenthesizedExpressionRule() {
		return getXParenthesizedExpressionAccess().getRule();
	}
	
	//XIfExpression XExpression:
	//	{XIfExpression}
	//	'if' '(' if=XExpression ')'
	//	then=XExpression (=> 'else' else=XExpression)?;
	public XbaseGrammarAccess.XIfExpressionElements getXIfExpressionAccess() {
		return gaXbase.getXIfExpressionAccess();
	}
	
	public ParserRule getXIfExpressionRule() {
		return getXIfExpressionAccess().getRule();
	}
	
	//XSwitchExpression XExpression:
	//	{XSwitchExpression}
	//	'switch' (=> ('(' declaredParam=JvmFormalParameter ':') switch=XExpression ')'
	//	| => (declaredParam=JvmFormalParameter ':')? switch=XExpression) '{'
	//	cases+=XCasePart* ('default' ':' default=XExpression)?
	//	'}';
	public XbaseGrammarAccess.XSwitchExpressionElements getXSwitchExpressionAccess() {
		return gaXbase.getXSwitchExpressionAccess();
	}
	
	public ParserRule getXSwitchExpressionRule() {
		return getXSwitchExpressionAccess().getRule();
	}
	
	//XCasePart:
	//	{XCasePart} typeGuard=JvmTypeReference? ('case' case=XExpression)? (':' then=XExpression | fallThrough?=',');
	public XbaseGrammarAccess.XCasePartElements getXCasePartAccess() {
		return gaXbase.getXCasePartAccess();
	}
	
	public ParserRule getXCasePartRule() {
		return getXCasePartAccess().getRule();
	}
	
	//XForLoopExpression XExpression:
	//	=> ({XForLoopExpression}
	//	'for' '(' declaredParam=JvmFormalParameter ':') forExpression=XExpression ')'
	//	eachExpression=XExpression;
	public XbaseGrammarAccess.XForLoopExpressionElements getXForLoopExpressionAccess() {
		return gaXbase.getXForLoopExpressionAccess();
	}
	
	public ParserRule getXForLoopExpressionRule() {
		return getXForLoopExpressionAccess().getRule();
	}
	
	//XBasicForLoopExpression XExpression:
	//	{XBasicForLoopExpression}
	//	'for' '(' (initExpressions+=XExpressionOrVarDeclaration (',' initExpressions+=XExpressionOrVarDeclaration)*)? ';'
	//	expression=XExpression? ';' (updateExpressions+=XExpression (',' updateExpressions+=XExpression)*)? ')'
	//	eachExpression=XExpression;
	public XbaseGrammarAccess.XBasicForLoopExpressionElements getXBasicForLoopExpressionAccess() {
		return gaXbase.getXBasicForLoopExpressionAccess();
	}
	
	public ParserRule getXBasicForLoopExpressionRule() {
		return getXBasicForLoopExpressionAccess().getRule();
	}
	
	//XWhileExpression XExpression:
	//	{XWhileExpression}
	//	'while' '(' predicate=XExpression ')'
	//	body=XExpression;
	public XbaseGrammarAccess.XWhileExpressionElements getXWhileExpressionAccess() {
		return gaXbase.getXWhileExpressionAccess();
	}
	
	public ParserRule getXWhileExpressionRule() {
		return getXWhileExpressionAccess().getRule();
	}
	
	//XDoWhileExpression XExpression:
	//	{XDoWhileExpression}
	//	'do'
	//	body=XExpression
	//	'while' '(' predicate=XExpression ')';
	public XbaseGrammarAccess.XDoWhileExpressionElements getXDoWhileExpressionAccess() {
		return gaXbase.getXDoWhileExpressionAccess();
	}
	
	public ParserRule getXDoWhileExpressionRule() {
		return getXDoWhileExpressionAccess().getRule();
	}
	
	//XBlockExpression XExpression:
	//	{XBlockExpression}
	//	'{' (expressions+=XExpressionOrVarDeclaration ';'?)*
	//	'}';
	public XbaseGrammarAccess.XBlockExpressionElements getXBlockExpressionAccess() {
		return gaXbase.getXBlockExpressionAccess();
	}
	
	public ParserRule getXBlockExpressionRule() {
		return getXBlockExpressionAccess().getRule();
	}
	
	//XExpressionOrVarDeclaration XExpression:
	//	XVariableDeclaration | XExpression;
	public XbaseGrammarAccess.XExpressionOrVarDeclarationElements getXExpressionOrVarDeclarationAccess() {
		return gaXbase.getXExpressionOrVarDeclarationAccess();
	}
	
	public ParserRule getXExpressionOrVarDeclarationRule() {
		return getXExpressionOrVarDeclarationAccess().getRule();
	}
	
	//XVariableDeclaration XExpression:
	//	{XVariableDeclaration} (writeable?='var' | 'val') (=> (type=JvmTypeReference name=ValidID) | name=ValidID) ('='
	//	right=XExpression)?;
	public XbaseGrammarAccess.XVariableDeclarationElements getXVariableDeclarationAccess() {
		return gaXbase.getXVariableDeclarationAccess();
	}
	
	public ParserRule getXVariableDeclarationRule() {
		return getXVariableDeclarationAccess().getRule();
	}
	
	//JvmFormalParameter types::JvmFormalParameter:
	//	parameterType=JvmTypeReference? name=ValidID;
	public XbaseGrammarAccess.JvmFormalParameterElements getJvmFormalParameterAccess() {
		return gaXbase.getJvmFormalParameterAccess();
	}
	
	public ParserRule getJvmFormalParameterRule() {
		return getJvmFormalParameterAccess().getRule();
	}
	
	//FullJvmFormalParameter types::JvmFormalParameter:
	//	parameterType=JvmTypeReference name=ValidID;
	public XbaseGrammarAccess.FullJvmFormalParameterElements getFullJvmFormalParameterAccess() {
		return gaXbase.getFullJvmFormalParameterAccess();
	}
	
	public ParserRule getFullJvmFormalParameterRule() {
		return getFullJvmFormalParameterAccess().getRule();
	}
	
	//XFeatureCall XExpression:
	//	{XFeatureCall} ('<' typeArguments+=JvmArgumentTypeReference (',' typeArguments+=JvmArgumentTypeReference)* '>')?
	//	feature=[types::JvmIdentifiableElement|IdOrSuper] (=> explicitOperationCall?='(' (featureCallArguments+=XShortClosure
	//	| featureCallArguments+=XExpression (',' featureCallArguments+=XExpression)*)?
	//	')')?
	//	featureCallArguments+=XClosure?;
	public XbaseGrammarAccess.XFeatureCallElements getXFeatureCallAccess() {
		return gaXbase.getXFeatureCallAccess();
	}
	
	public ParserRule getXFeatureCallRule() {
		return getXFeatureCallAccess().getRule();
	}
	
	//FeatureCallID:
	//	ValidID | 'extends' | 'static' | 'import' | 'extension';
	public XbaseGrammarAccess.FeatureCallIDElements getFeatureCallIDAccess() {
		return gaXbase.getFeatureCallIDAccess();
	}
	
	public ParserRule getFeatureCallIDRule() {
		return getFeatureCallIDAccess().getRule();
	}
	
	//IdOrSuper:
	//	FeatureCallID | 'super';
	public XbaseGrammarAccess.IdOrSuperElements getIdOrSuperAccess() {
		return gaXbase.getIdOrSuperAccess();
	}
	
	public ParserRule getIdOrSuperRule() {
		return getIdOrSuperAccess().getRule();
	}
	
	//XConstructorCall XExpression:
	//	{XConstructorCall}
	//	'new' constructor=[types::JvmConstructor|QualifiedName] (=> '<' typeArguments+=JvmArgumentTypeReference (','
	//	typeArguments+=JvmArgumentTypeReference)* '>')? (=> explicitConstructorCall?='(' (arguments+=XShortClosure
	//	| arguments+=XExpression (',' arguments+=XExpression)*)?
	//	')')?
	//	arguments+=XClosure?;
	public XbaseGrammarAccess.XConstructorCallElements getXConstructorCallAccess() {
		return gaXbase.getXConstructorCallAccess();
	}
	
	public ParserRule getXConstructorCallRule() {
		return getXConstructorCallAccess().getRule();
	}
	
	//XBooleanLiteral XExpression:
	//	{XBooleanLiteral} ('false' | isTrue?='true');
	public XbaseGrammarAccess.XBooleanLiteralElements getXBooleanLiteralAccess() {
		return gaXbase.getXBooleanLiteralAccess();
	}
	
	public ParserRule getXBooleanLiteralRule() {
		return getXBooleanLiteralAccess().getRule();
	}
	
	//XNullLiteral XExpression:
	//	{XNullLiteral} 'null';
	public XbaseGrammarAccess.XNullLiteralElements getXNullLiteralAccess() {
		return gaXbase.getXNullLiteralAccess();
	}
	
	public ParserRule getXNullLiteralRule() {
		return getXNullLiteralAccess().getRule();
	}
	
	//XNumberLiteral XExpression:
	//	{XNumberLiteral} value=Number;
	public XbaseGrammarAccess.XNumberLiteralElements getXNumberLiteralAccess() {
		return gaXbase.getXNumberLiteralAccess();
	}
	
	public ParserRule getXNumberLiteralRule() {
		return getXNumberLiteralAccess().getRule();
	}
	
	//XStringLiteral XExpression:
	//	{XStringLiteral} value=STRING;
	public XbaseGrammarAccess.XStringLiteralElements getXStringLiteralAccess() {
		return gaXbase.getXStringLiteralAccess();
	}
	
	public ParserRule getXStringLiteralRule() {
		return getXStringLiteralAccess().getRule();
	}
	
	//XTypeLiteral XExpression:
	//	{XTypeLiteral} 'typeof' '(' type=[types::JvmType|QualifiedName] arrayDimensions+=ArrayBrackets* ')';
	public XbaseGrammarAccess.XTypeLiteralElements getXTypeLiteralAccess() {
		return gaXbase.getXTypeLiteralAccess();
	}
	
	public ParserRule getXTypeLiteralRule() {
		return getXTypeLiteralAccess().getRule();
	}
	
	//XThrowExpression XExpression:
	//	{XThrowExpression} 'throw' expression=XExpression;
	public XbaseGrammarAccess.XThrowExpressionElements getXThrowExpressionAccess() {
		return gaXbase.getXThrowExpressionAccess();
	}
	
	public ParserRule getXThrowExpressionRule() {
		return getXThrowExpressionAccess().getRule();
	}
	
	//XReturnExpression XExpression:
	//	{XReturnExpression} 'return' -> expression=XExpression?;
	public XbaseGrammarAccess.XReturnExpressionElements getXReturnExpressionAccess() {
		return gaXbase.getXReturnExpressionAccess();
	}
	
	public ParserRule getXReturnExpressionRule() {
		return getXReturnExpressionAccess().getRule();
	}
	
	//XTryCatchFinallyExpression XExpression:
	//	{XTryCatchFinallyExpression}
	//	'try'
	//	expression=XExpression (catchClauses+=XCatchClause+ (=> 'finally' finallyExpression=XExpression)?
	//	| 'finally' finallyExpression=XExpression);
	public XbaseGrammarAccess.XTryCatchFinallyExpressionElements getXTryCatchFinallyExpressionAccess() {
		return gaXbase.getXTryCatchFinallyExpressionAccess();
	}
	
	public ParserRule getXTryCatchFinallyExpressionRule() {
		return getXTryCatchFinallyExpressionAccess().getRule();
	}
	
	//XSynchronizedExpression XExpression:
	//	=> ({XSynchronizedExpression}
	//	'synchronized' '(') param=XExpression ')' expression=XExpression;
	public XbaseGrammarAccess.XSynchronizedExpressionElements getXSynchronizedExpressionAccess() {
		return gaXbase.getXSynchronizedExpressionAccess();
	}
	
	public ParserRule getXSynchronizedExpressionRule() {
		return getXSynchronizedExpressionAccess().getRule();
	}
	
	//XCatchClause:
	//	=> 'catch' '(' declaredParam=FullJvmFormalParameter ')' expression=XExpression;
	public XbaseGrammarAccess.XCatchClauseElements getXCatchClauseAccess() {
		return gaXbase.getXCatchClauseAccess();
	}
	
	public ParserRule getXCatchClauseRule() {
		return getXCatchClauseAccess().getRule();
	}
	
	//@Override
	//QualifiedName:
	//	ValidID (=> '.' ValidID)*;
	public XbaseGrammarAccess.QualifiedNameElements getQualifiedNameAccess() {
		return gaXbase.getQualifiedNameAccess();
	}
	
	public ParserRule getQualifiedNameRule() {
		return getQualifiedNameAccess().getRule();
	}
	
	//Number hidden():
	//	HEX | (INT | DECIMAL) ('.' (INT | DECIMAL))?;
	public XbaseGrammarAccess.NumberElements getNumberAccess() {
		return gaXbase.getNumberAccess();
	}
	
	public ParserRule getNumberRule() {
		return getNumberAccess().getRule();
	}
	
	///**
	// * Dummy rule, for "better" downwards compatibility, since GrammarAccess generates non-static inner classes,
	// * which makes downstream grammars break on classloading, when a rule is removed.
	// */ StaticQualifier:
	//	(ValidID '::')+;
	public XbaseGrammarAccess.StaticQualifierElements getStaticQualifierAccess() {
		return gaXbase.getStaticQualifierAccess();
	}
	
	public ParserRule getStaticQualifierRule() {
		return getStaticQualifierAccess().getRule();
	}
	
	//terminal HEX:
	//	('0x' | '0X') ('0'..'9' | 'a'..'f' | 'A'..'F' | '_')+ ('#' (('b' | 'B') ('i' | 'I') | ('l' | 'L')))?;
	public TerminalRule getHEXRule() {
		return gaXbase.getHEXRule();
	}
	
	//terminal INT returns ecore::EInt:
	//	'0'..'9' ('0'..'9' | '_')*;
	public TerminalRule getINTRule() {
		return gaXbase.getINTRule();
	}
	
	//terminal DECIMAL:
	//	INT (('e' | 'E') ('+' | '-')? INT)? (('b' | 'B') ('i' | 'I' | 'd' | 'D') | ('l' | 'L' | 'd' | 'D' | 'f' | 'F'))?;
	public TerminalRule getDECIMALRule() {
		return gaXbase.getDECIMALRule();
	}
	
	//JvmTypeReference:
	//	JvmParameterizedTypeReference => ({JvmGenericArrayTypeReference.componentType=current} ArrayBrackets)*
	//	| XFunctionTypeRef;
	public XtypeGrammarAccess.JvmTypeReferenceElements getJvmTypeReferenceAccess() {
		return gaXtype.getJvmTypeReferenceAccess();
	}
	
	public ParserRule getJvmTypeReferenceRule() {
		return getJvmTypeReferenceAccess().getRule();
	}
	
	//ArrayBrackets:
	//	'[' ']';
	public XtypeGrammarAccess.ArrayBracketsElements getArrayBracketsAccess() {
		return gaXtype.getArrayBracketsAccess();
	}
	
	public ParserRule getArrayBracketsRule() {
		return getArrayBracketsAccess().getRule();
	}
	
	//XFunctionTypeRef:
	//	('(' (paramTypes+=JvmTypeReference (',' paramTypes+=JvmTypeReference)*)? ')')? '=>' returnType=JvmTypeReference;
	public XtypeGrammarAccess.XFunctionTypeRefElements getXFunctionTypeRefAccess() {
		return gaXtype.getXFunctionTypeRefAccess();
	}
	
	public ParserRule getXFunctionTypeRefRule() {
		return getXFunctionTypeRefAccess().getRule();
	}
	
	//JvmParameterizedTypeReference:
	//	type=[JvmType|super::QualifiedName] (=> '<' arguments+=JvmArgumentTypeReference (','
	//	arguments+=JvmArgumentTypeReference)* '>' (=> ({JvmInnerTypeReference.outer=current} '.') type=[JvmType|ValidID] (=>
	//	'<' arguments+=JvmArgumentTypeReference (',' arguments+=JvmArgumentTypeReference)* '>')?)*)?;
	public XtypeGrammarAccess.JvmParameterizedTypeReferenceElements getJvmParameterizedTypeReferenceAccess() {
		return gaXtype.getJvmParameterizedTypeReferenceAccess();
	}
	
	public ParserRule getJvmParameterizedTypeReferenceRule() {
		return getJvmParameterizedTypeReferenceAccess().getRule();
	}
	
	//JvmArgumentTypeReference JvmTypeReference:
	//	JvmTypeReference | JvmWildcardTypeReference;
	public XtypeGrammarAccess.JvmArgumentTypeReferenceElements getJvmArgumentTypeReferenceAccess() {
		return gaXtype.getJvmArgumentTypeReferenceAccess();
	}
	
	public ParserRule getJvmArgumentTypeReferenceRule() {
		return getJvmArgumentTypeReferenceAccess().getRule();
	}
	
	//JvmWildcardTypeReference:
	//	{JvmWildcardTypeReference} '?' (constraints+=JvmUpperBound constraints+=JvmUpperBoundAnded*
	//	| constraints+=JvmLowerBound constraints+=JvmLowerBoundAnded*)?;
	public XtypeGrammarAccess.JvmWildcardTypeReferenceElements getJvmWildcardTypeReferenceAccess() {
		return gaXtype.getJvmWildcardTypeReferenceAccess();
	}
	
	public ParserRule getJvmWildcardTypeReferenceRule() {
		return getJvmWildcardTypeReferenceAccess().getRule();
	}
	
	//JvmUpperBound:
	//	'extends' typeReference=JvmTypeReference;
	public XtypeGrammarAccess.JvmUpperBoundElements getJvmUpperBoundAccess() {
		return gaXtype.getJvmUpperBoundAccess();
	}
	
	public ParserRule getJvmUpperBoundRule() {
		return getJvmUpperBoundAccess().getRule();
	}
	
	//JvmUpperBoundAnded JvmUpperBound:
	//	'&' typeReference=JvmTypeReference;
	public XtypeGrammarAccess.JvmUpperBoundAndedElements getJvmUpperBoundAndedAccess() {
		return gaXtype.getJvmUpperBoundAndedAccess();
	}
	
	public ParserRule getJvmUpperBoundAndedRule() {
		return getJvmUpperBoundAndedAccess().getRule();
	}
	
	//JvmLowerBound:
	//	'super' typeReference=JvmTypeReference;
	public XtypeGrammarAccess.JvmLowerBoundElements getJvmLowerBoundAccess() {
		return gaXtype.getJvmLowerBoundAccess();
	}
	
	public ParserRule getJvmLowerBoundRule() {
		return getJvmLowerBoundAccess().getRule();
	}
	
	//JvmLowerBoundAnded JvmLowerBound:
	//	'&' typeReference=JvmTypeReference;
	public XtypeGrammarAccess.JvmLowerBoundAndedElements getJvmLowerBoundAndedAccess() {
		return gaXtype.getJvmLowerBoundAndedAccess();
	}
	
	public ParserRule getJvmLowerBoundAndedRule() {
		return getJvmLowerBoundAndedAccess().getRule();
	}
	
	//JvmTypeParameter:
	//	name=ValidID (constraints+=JvmUpperBound constraints+=JvmUpperBoundAnded*)?;
	public XtypeGrammarAccess.JvmTypeParameterElements getJvmTypeParameterAccess() {
		return gaXtype.getJvmTypeParameterAccess();
	}
	
	public ParserRule getJvmTypeParameterRule() {
		return getJvmTypeParameterAccess().getRule();
	}
	
	//QualifiedNameWithWildcard:
	//	super::QualifiedName '.' '*';
	public XtypeGrammarAccess.QualifiedNameWithWildcardElements getQualifiedNameWithWildcardAccess() {
		return gaXtype.getQualifiedNameWithWildcardAccess();
	}
	
	public ParserRule getQualifiedNameWithWildcardRule() {
		return getQualifiedNameWithWildcardAccess().getRule();
	}
	
	//ValidID:
	//	ID;
	public XtypeGrammarAccess.ValidIDElements getValidIDAccess() {
		return gaXtype.getValidIDAccess();
	}
	
	public ParserRule getValidIDRule() {
		return getValidIDAccess().getRule();
	}
	
	//XImportSection:
	//	importDeclarations+=XImportDeclaration+;
	public XtypeGrammarAccess.XImportSectionElements getXImportSectionAccess() {
		return gaXtype.getXImportSectionAccess();
	}
	
	public ParserRule getXImportSectionRule() {
		return getXImportSectionAccess().getRule();
	}
	
	//XImportDeclaration:
	//	'import' (static?='static' extension?='extension'? importedType=[JvmDeclaredType|QualifiedNameInStaticImport]
	//	(wildcard?='*' | memberName=ValidID) | importedType=[JvmDeclaredType|super::QualifiedName] |
	//	importedNamespace=QualifiedNameWithWildcard) ';'?;
	public XtypeGrammarAccess.XImportDeclarationElements getXImportDeclarationAccess() {
		return gaXtype.getXImportDeclarationAccess();
	}
	
	public ParserRule getXImportDeclarationRule() {
		return getXImportDeclarationAccess().getRule();
	}
	
	//QualifiedNameInStaticImport:
	//	(ValidID '.')+;
	public XtypeGrammarAccess.QualifiedNameInStaticImportElements getQualifiedNameInStaticImportAccess() {
		return gaXtype.getQualifiedNameInStaticImportAccess();
	}
	
	public ParserRule getQualifiedNameInStaticImportRule() {
		return getQualifiedNameInStaticImportAccess().getRule();
	}
	
	//terminal ID:
	//	'^'? ('a'..'z' | 'A'..'Z' | '$' | '_') ('a'..'z' | 'A'..'Z' | '$' | '_' | '0'..'9')*;
	public TerminalRule getIDRule() {
		return gaXtype.getIDRule();
	}
	
	//terminal STRING:
	//	'"' ('\\' . | !('\\' | '"'))* '"'? |
	//	"'" ('\\' . | !('\\' | "'"))* "'"?;
	public TerminalRule getSTRINGRule() {
		return gaXtype.getSTRINGRule();
	}
	
	//terminal ML_COMMENT:
	//	'/*'->'*/';
	public TerminalRule getML_COMMENTRule() {
		return gaXtype.getML_COMMENTRule();
	}
	
	//terminal SL_COMMENT:
	//	'//' !('\n' | '\r')* ('\r'? '\n')?;
	public TerminalRule getSL_COMMENTRule() {
		return gaXtype.getSL_COMMENTRule();
	}
	
	//terminal WS:
	//	' ' | '\t' | '\r' | '\n'+;
	public TerminalRule getWSRule() {
		return gaXtype.getWSRule();
	}
	
	//terminal ANY_OTHER:
	//	.;
	public TerminalRule getANY_OTHERRule() {
		return gaXtype.getANY_OTHERRule();
	}
}
