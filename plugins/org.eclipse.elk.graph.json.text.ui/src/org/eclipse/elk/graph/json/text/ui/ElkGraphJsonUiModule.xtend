/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 *
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.graph.json.text.ui

import org.eclipse.xtend.lib.annotations.FinalFieldsConstructor

/**
 * Use this class to register components to be used within the Eclipse IDE.
 */
@FinalFieldsConstructor
class ElkGraphJsonUiModule extends AbstractElkGraphJsonUiModule {
    public static val String PLUGIN_ID = "org.eclipse.elk.graph.json.text.ui";
}
