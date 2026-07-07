/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.meta.ui.labeling

import com.google.inject.Inject
import org.eclipse.elk.core.meta.metaData.MdAlgorithm
import org.eclipse.elk.core.meta.metaData.MdCategory
import org.eclipse.elk.core.meta.metaData.MdGroup
import org.eclipse.elk.core.meta.metaData.MdOption
import org.eclipse.emf.edit.ui.provider.AdapterFactoryLabelProvider
import org.eclipse.xtext.ui.label.DefaultEObjectLabelProvider

/**
 * Provides labels for EObjects.
 * 
 * See https://www.eclipse.org/Xtext/documentation/304_ide_concepts.html#label-provider
 */
class MetaDataLabelProvider extends DefaultEObjectLabelProvider {

	@Inject
	new(AdapterFactoryLabelProvider delegate) {
		super(delegate);
	}

	
	def text(MdAlgorithm algorithm) {
	    'Algorithm: ' + algorithm.name;
	}
	
	def text(MdCategory category) {
	    'Category: ' + category.name;
	}
	
	def text(MdGroup group) {
	    'Group: ' + group.name;
	}
	
	def text(MdOption option) {
	    'Option: ' + option.name;
	}
	
	// Labels and icons can be computed like this:
//
//	def image(Greeting ele) {
//		'Greeting.gif'
//	}
}
