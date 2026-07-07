/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 *
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.graph.json.text


/**
 * Initialization support for running Xtext languages without Equinox extension registry.
 */
class ElkGraphJsonStandaloneSetup extends ElkGraphJsonStandaloneSetupGenerated {

	def static void doSetup() {
		new ElkGraphJsonStandaloneSetup().createInjectorAndDoEMFRegistration()
	}
}
