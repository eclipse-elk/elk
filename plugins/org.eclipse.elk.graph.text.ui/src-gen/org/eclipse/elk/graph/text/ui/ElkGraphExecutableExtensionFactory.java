/*******************************************************************************
 * Copyright (c) 2016 TypeFox GmbH (http://www.typefox.io) and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.graph.text.ui;

import com.google.inject.Injector;
import org.eclipse.elk.graph.text.ui.internal.TextActivator;
import org.eclipse.xtext.ui.guice.AbstractGuiceAwareExecutableExtensionFactory;
import org.osgi.framework.Bundle;

/**
 * This class was generated. Customizations should only happen in a newly
 * introduced subclass. 
 */
public class ElkGraphExecutableExtensionFactory extends AbstractGuiceAwareExecutableExtensionFactory {

	@Override
	protected Bundle getBundle() {
		return TextActivator.getInstance().getBundle();
	}
	
	@Override
	protected Injector getInjector() {
		return TextActivator.getInstance().getInjector(TextActivator.ORG_ECLIPSE_ELK_GRAPH_TEXT_ELKGRAPH);
	}
	
}
