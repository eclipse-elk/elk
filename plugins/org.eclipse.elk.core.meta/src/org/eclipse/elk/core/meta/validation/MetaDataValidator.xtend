/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.meta.validation

import java.util.Map
import org.eclipse.elk.core.meta.metaData.MdAlgorithm
import org.eclipse.elk.core.meta.metaData.MdBundle
import org.eclipse.elk.core.meta.metaData.MdBundleMember
import org.eclipse.elk.core.meta.metaData.MdCategory
import org.eclipse.elk.core.meta.metaData.MdGroup
import org.eclipse.elk.core.meta.metaData.MdOption
import org.eclipse.xtext.validation.Check

import static org.eclipse.elk.core.meta.metaData.MetaDataPackage.Literals.*

/**
 * This class contains custom validation rules. 
 *
 * See https://www.eclipse.org/Xtext/documentation/303_runtime_concepts.html#validation
 */
class MetaDataValidator extends AbstractMetaDataValidator {
	
    @Check
    def void checkDuplicateMemberId(MdBundle bundle) {
        bundle.members.checkDuplicateIds
    }
    
    private def void checkDuplicateIds(Iterable<? extends MdBundleMember> elements) {
        val Map<String, MdAlgorithm> algorithmIds = newHashMap
        val Map<String, MdCategory> categoryIds = newHashMap
        val Map<String, MdOption> propertyIds = newHashMap
        val Map<String, MdGroup> groupIds = newHashMap
        for (element : elements) {
            switch element {
                MdAlgorithm: algorithmIds.checkExistsAndRemember(element)
                MdCategory: categoryIds.checkExistsAndRemember(element)
                MdGroup: {
                    groupIds.checkExistsAndRemember(element)
                    element.children.checkDuplicateIds
                }
                MdOption: propertyIds.checkExistsAndRemember(element)
            }
        }
    }
    
    private def  <T extends MdBundleMember> void checkExistsAndRemember(Map<String, T> map, T element) {
        if (map.containsKey(element.name)) {
                val otherMember = map.get(element.name)
                if (otherMember !== null) {
                    duplicateId(otherMember)
                    // The first occurrence should be marked only once
                    map.put(element.name, null)
                }
                duplicateId(element)
            } else {
                map.put(element.name, element)
            }
    }
    
    private def void duplicateId(MdBundleMember member) {
        error("The id '" + member.name + "' is already used.", member, MD_BUNDLE_MEMBER__NAME)
    }
	
}
