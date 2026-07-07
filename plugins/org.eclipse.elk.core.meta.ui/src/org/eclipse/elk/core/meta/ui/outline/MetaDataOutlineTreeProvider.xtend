/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.meta.ui.outline

import org.eclipse.elk.core.meta.metaData.MdAlgorithm
import org.eclipse.elk.core.meta.metaData.MdBundle
import org.eclipse.elk.core.meta.metaData.MdBundleMember
import org.eclipse.elk.core.meta.metaData.MdCategory
import org.eclipse.elk.core.meta.metaData.MdGroup
import org.eclipse.elk.core.meta.metaData.MdGroupOrOption
import org.eclipse.elk.core.meta.metaData.MdModel
import org.eclipse.elk.core.meta.metaData.MdOption
import org.eclipse.xtext.ui.editor.outline.IOutlineNode
import org.eclipse.xtext.ui.editor.outline.impl.DefaultOutlineTreeProvider
import org.eclipse.xtext.ui.editor.outline.impl.DocumentRootNode

/**
 * Customization of the default outline structure.
 *
 * See https://www.eclipse.org/Xtext/documentation/304_ide_concepts.html#outline
 */
class MetaDataOutlineTreeProvider extends DefaultOutlineTreeProvider {
    
    protected def void _createChildren(DocumentRootNode rootNode, MdModel model) {
        if (model.bundle !== null) {
            createNode(rootNode, model.bundle);
        }
    }
    
    protected def void _createChildren(IOutlineNode parentNode, MdBundle bundle) {
        for (MdBundleMember member : bundle.members) {
            if (member instanceof MdCategory
                || member instanceof MdAlgorithm
                || member instanceof MdOption
                || member instanceof MdGroup) {
                    
                createNode(parentNode, member);
            }
        }
    }
    
    protected def void _createChildren(IOutlineNode parentNode, MdCategory category) {
        // Categories shouldn't have children
        return;
    }
    
    protected def boolean _isLeaf(MdCategory category) {
        return true;
    }
    
    protected def void _createChildren(IOutlineNode parentNode, MdAlgorithm algorithm) {
        // Algorithms shouldn't have children
        return;
    }
    
    protected def boolean _isLeaf(MdAlgorithm algorithm) {
        return true;
    }
    
    protected def void _createChildren(IOutlineNode parentNode, MdOption option) {
        // Options shouldn't have children
        return;
    }
    
    protected def boolean _isLeaf(MdOption option) {
        return true;
    }
    
    protected def void _createChildren(IOutlineNode parentNode, MdGroup group) {
        for (MdGroupOrOption child : group.children) {
            createNode(parentNode, child);
        }
    }
    
}
