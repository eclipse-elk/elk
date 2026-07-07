/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 *
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.graph.json.text.ide

import org.eclipse.elk.graph.json.text.ide.contentassist.ElkGraphJsonProposalProvider
import org.eclipse.xtext.ide.editor.contentassist.IdeContentProposalProvider

/**
 * Use this class to register ide components.
 */
class ElkGraphJsonIdeModule extends AbstractElkGraphJsonIdeModule {
        
    def Class<? extends IdeContentProposalProvider> bindIdeContentProposalProvider() {
        ElkGraphJsonProposalProvider
    }
    
}
