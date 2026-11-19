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
package org.eclipse.elk.core.meta.formatting2

import org.eclipse.xtext.formatting2.AbstractFormatter2
import org.eclipse.xtext.formatting2.IFormattableDocument
import org.eclipse.elk.core.meta.metaData.MdBundle
import org.eclipse.xtext.xbase.lib.Procedures.Procedure1
import org.eclipse.xtext.formatting2.IHiddenRegionFormatter
import org.eclipse.elk.core.meta.services.MetaDataGrammarAccess
import com.google.inject.Inject
import org.eclipse.elk.core.meta.metaData.MdGroup
import org.eclipse.elk.core.meta.metaData.MdOption
import org.eclipse.elk.core.meta.metaData.MdOptionDependency
import org.eclipse.elk.core.meta.metaData.MdAlgorithm
import org.eclipse.elk.core.meta.metaData.MdCategory
import org.eclipse.elk.core.meta.metaData.MdOptionSupport
import org.eclipse.elk.core.meta.metaData.MdBundleMember
import org.eclipse.elk.core.meta.metaData.MdGroupOrOption

class MetaDataFormatter extends AbstractFormatter2 {

    static val Procedure1<? super IHiddenRegionFormatter> no_space = [noSpace]
    static val Procedure1<? super IHiddenRegionFormatter> one_space = [oneSpace]
    static val Procedure1<? super IHiddenRegionFormatter> new_line = [newLine]
    static val Procedure1<? super IHiddenRegionFormatter> indention = [indent]
    static val Procedure1<? super IHiddenRegionFormatter> new_line_twice = [setNewLines(2, 2, 2)]

    @Inject extension MetaDataGrammarAccess

    def dispatch void format(MdBundleMember bundleMember, extension IFormattableDocument document) {
        if (bundleMember instanceof MdAlgorithm) {
            bundleMember.format
        } else if (bundleMember instanceof MdCategory) {
            bundleMember.format
        } else if (bundleMember instanceof MdGroupOrOption) {
            bundleMember.format
        }
    }

    def dispatch void format(MdGroupOrOption grOpt, extension IFormattableDocument document) {
        if (grOpt instanceof MdGroup) {
            grOpt.format
        } else if (grOpt instanceof MdOption) {
            grOpt.format
        }
    }

    def dispatch void format(MdBundle bundle, extension IFormattableDocument document) {
        bundle.regionFor.keyword(mdBundleAccess.bundleKeyword_1_0).append(one_space)
        interior(
            bundle.regionFor.keyword(mdBundleAccess.leftCurlyBracketKeyword_1_1).prepend(one_space),
            bundle.regionFor.keyword(mdBundleAccess.rightCurlyBracketKeyword_1_3),
            indention
        )

        bundle.regionFor.keyword(mdBundleAccess.labelKeyword_1_2_0_0).prepend(new_line).append(one_space)
        bundle.regionFor.keyword(mdBundleAccess.metadataClassKeyword_1_2_1_0).prepend(new_line).append(one_space)
        bundle.regionFor.keyword(mdBundleAccess.documentationFolderKeyword_1_2_2_0).prepend(new_line).append(one_space)
        bundle.regionFor.keyword(mdBundleAccess.idPrefixKeyword_1_2_3_0).prepend(new_line).append(one_space)

        for (member : bundle.members) {
            member.prepend(new_line_twice).format
        }
    }

    def dispatch void format(MdGroup group, extension IFormattableDocument document) {
        group.regionFor.keyword(mdGroupAccess.groupKeyword_0).append(one_space)
        interior(
            group.regionFor.keyword(mdGroupAccess.leftCurlyBracketKeyword_2).prepend(one_space),
            group.regionFor.keyword(mdGroupAccess.rightCurlyBracketKeyword_5).prepend(new_line_twice),
            indention
        )

        group.regionFor.keyword(mdGroupAccess.documentationKeyword_3_0).append(one_space)

        for (child : group.children) {
            child.prepend(new_line_twice).format
        }
    }

    def dispatch void format(MdOption option, extension IFormattableDocument document) {
        option.regionFor.keyword(mdOptionAccess.deprecatedDeprecatedKeyword_0_0).append(one_space)
        option.regionFor.keyword(mdOptionAccess.advancedAdvancedKeyword_1_0_0).append(one_space)
        option.regionFor.keyword(mdOptionAccess.programmaticProgrammaticKeyword_1_1_0).append(one_space)
        option.regionFor.keyword(mdOptionAccess.outputOutputKeyword_1_2_0).append(one_space)
        option.regionFor.keyword(mdOptionAccess.globalGlobalKeyword_1_3_0).append(one_space)

        option.regionFor.keyword(mdOptionAccess.optionKeyword_2).append(one_space)
        interior(
            option.regionFor.keyword(mdOptionAccess.leftCurlyBracketKeyword_5).prepend(one_space).append(new_line),
            option.regionFor.keyword(mdOptionAccess.rightCurlyBracketKeyword_8),
            indention
        ) 

        option.regionFor.keyword(mdOptionAccess.labelKeyword_6_0_0).prepend(new_line).append(one_space)
        option.regionFor.keyword(mdOptionAccess.descriptionKeyword_6_1_0).prepend(new_line).append(new_line)
        interior(
            option.regionFor.keyword(mdOptionAccess.descriptionKeyword_6_1_0),
            option.regionFor.assignment(mdOptionAccess.descriptionAssignment_6_1_1).nextSemanticRegion,
            indention
        )
        option.regionFor.keyword(mdOptionAccess.documentationKeyword_6_2_0).prepend(new_line).append(one_space)
        option.regionFor.keyword(mdOptionAccess.defaultKeyword_6_3_0).prepend(new_line).append(one_space)
        option.regionFor.keyword(mdOptionAccess.equalsSignKeyword_6_3_1).append(one_space)
        option.regionFor.keyword(mdOptionAccess.lowerBoundKeyword_6_4_0).prepend(new_line).append(one_space)
        option.regionFor.keyword(mdOptionAccess.equalsSignKeyword_6_4_1).append(one_space)
        option.regionFor.keyword(mdOptionAccess.upperBoundKeyword_6_5_0).prepend(new_line).append(one_space)
        option.regionFor.keyword(mdOptionAccess.equalsSignKeyword_6_5_1).append(one_space)
        option.regionFor.keyword(mdOptionAccess.targetsKeyword_6_6_0).prepend(new_line).append(one_space)
        option.regionFor.keyword(mdOptionAccess.commaKeyword_6_6_2_0).prepend(no_space).append(one_space)
        option.regionFor.keyword(mdOptionAccess.legacyIdsKeyword_6_7_0).prepend(new_line).append(one_space)
        option.regionFor.keyword(mdOptionAccess.commaKeyword_6_7_2_0).prepend(no_space).append(one_space)

        for (dependencies : option.dependencies) {
            dependencies.append(new_line)
        }
    }

    def dispatch void format(MdOptionDependency optionDependency, extension IFormattableDocument document) {
        optionDependency.regionFor.keyword(mdOptionDependencyAccess.requiresKeyword_0).append(one_space)
        optionDependency.regionFor.keyword(mdOptionDependencyAccess.equalsSignEqualsSignKeyword_2_0).prepend(one_space).
            append(one_space)
    }

    def dispatch void format(MdAlgorithm algorithm, extension IFormattableDocument document) {
        algorithm.regionFor.keyword(mdAlgorithmAccess.deprecatedDeprecatedKeyword_0_0).append(one_space)

        algorithm.regionFor.keyword(mdAlgorithmAccess.algorithmKeyword_1).append(one_space)
        algorithm.regionFor.keyword(mdAlgorithmAccess.leftParenthesisKeyword_3).prepend(no_space)
        interior(
            algorithm.regionFor.keyword(mdAlgorithmAccess.leftCurlyBracketKeyword_7).prepend(one_space).append(
                new_line),
            algorithm.regionFor.keyword(mdAlgorithmAccess.rightCurlyBracketKeyword_10),
            indention
        )

        algorithm.regionFor.keyword(mdAlgorithmAccess.labelKeyword_8_0_0).prepend(new_line).append(one_space)
        algorithm.regionFor.keyword(mdAlgorithmAccess.metadataClassKeyword_8_1_0).prepend(new_line).append(one_space)
        algorithm.regionFor.keyword(mdAlgorithmAccess.descriptionKeyword_8_2_0).prepend(new_line).append(new_line)
        interior(
            algorithm.regionFor.keyword(mdAlgorithmAccess.descriptionKeyword_8_2_0),
            algorithm.regionFor.assignment(mdAlgorithmAccess.descriptionAssignment_8_2_1).nextSemanticRegion,
            indention
        )
        algorithm.regionFor.keyword(mdAlgorithmAccess.documentationKeyword_8_3_0).prepend(new_line).append(one_space)
        algorithm.regionFor.keyword(mdAlgorithmAccess.categoryKeyword_8_4_0).prepend(new_line).append(one_space)
        algorithm.regionFor.keyword(mdAlgorithmAccess.previewKeyword_8_5_0).prepend(new_line).append(one_space)
        algorithm.regionFor.keyword(mdAlgorithmAccess.featuresKeyword_8_6_0).prepend(new_line).append(one_space)
        algorithm.regionFor.keyword(mdAlgorithmAccess.validatorKeyword_8_7_0).prepend(one_space).append(one_space)

        for (option : algorithm.supportedOptions) {
            option.append(new_line).format
        }
    }

    def dispatch void format(MdCategory category, extension IFormattableDocument document) {
        category.regionFor.keyword(mdCategoryAccess.deprecatedDeprecatedKeyword_0_0).append(one_space)

        category.regionFor.keyword(mdCategoryAccess.categoryKeyword_1).append(one_space)
        interior(
            category.regionFor.keyword(mdCategoryAccess.leftCurlyBracketKeyword_3).prepend(one_space).append(new_line),
            category.regionFor.keyword(mdCategoryAccess.rightCurlyBracketKeyword_5),
            indention
        )

        category.regionFor.keyword(mdCategoryAccess.labelKeyword_4_0_0).prepend(new_line).append(one_space)
        category.regionFor.keyword(mdCategoryAccess.descriptionKeyword_4_1_0).prepend(new_line).append(new_line)
        interior(
            category.regionFor.keyword(mdCategoryAccess.descriptionKeyword_4_1_0),
            category.regionFor.assignment(mdCategoryAccess.descriptionAssignment_4_1_1).nextSemanticRegion,
            indention
        )
        category.regionFor.keyword(mdCategoryAccess.documentationKeyword_4_2_0).prepend(new_line).append(one_space)
    }

    def dispatch void format(MdOptionSupport optionSupport, extension IFormattableDocument document) {
        optionSupport.regionFor.keyword(mdOptionSupportAccess.supportsKeyword_0).append(one_space)
        optionSupport.regionFor.keyword(mdOptionSupportAccess.equalsSignKeyword_2_0).prepend(one_space).append(
            one_space)
        optionSupport.regionFor.keyword(mdOptionSupportAccess.documentationKeyword_3_0).prepend(new_line).append(
            one_space)
    }
}
