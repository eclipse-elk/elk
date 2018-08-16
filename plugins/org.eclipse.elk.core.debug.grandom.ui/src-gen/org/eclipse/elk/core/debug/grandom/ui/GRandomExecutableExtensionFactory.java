/*******************************************************************************
 * Copyright (c) 2018 Kiel University and others.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 *******************************************************************************/
package org.eclipse.elk.core.debug.grandom.ui;

import com.google.inject.Injector;
import org.eclipse.elk.core.debug.grandom.ui.internal.GrandomActivator;
import org.eclipse.xtext.ui.guice.AbstractGuiceAwareExecutableExtensionFactory;
import org.osgi.framework.Bundle;

/**
 * This class was generated. Customizations should only happen in a newly
 * introduced subclass. 
 */
public class GRandomExecutableExtensionFactory extends AbstractGuiceAwareExecutableExtensionFactory {

	@Override
	protected Bundle getBundle() {
		return GrandomActivator.getInstance().getBundle();
	}
	
	@Override
	protected Injector getInjector() {
		return GrandomActivator.getInstance().getInjector(GrandomActivator.ORG_ECLIPSE_ELK_CORE_DEBUG_GRANDOM_GRANDOM);
	}
	
}
