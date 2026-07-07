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
package org.eclipse.elk.core.debug.grandom.ide

import com.google.inject.Guice
import org.eclipse.elk.core.debug.grandom.GRandomRuntimeModule
import org.eclipse.elk.core.debug.grandom.GRandomStandaloneSetup
import org.eclipse.xtext.util.Modules2

/**
 * Initialization support for running Xtext languages as language servers.
 */
class GRandomIdeSetup extends GRandomStandaloneSetup {

	override createInjector() {
		Guice.createInjector(Modules2.mixin(new GRandomRuntimeModule, new GRandomIdeModule))
	}
	
}
