/*******************************************************************************
 * Copyright (c) 2008, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.

 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.graphviz.layouter.preferences;

import org.eclipse.core.runtime.preferences.AbstractPreferenceInitializer;
import org.eclipse.elk.alg.graphviz.layouter.GraphvizLayoutProvider;
import org.eclipse.elk.alg.graphviz.layouter.GraphvizTool;
import org.eclipse.jface.preference.IPreferenceStore;

/**
 * Class used to initialize default preference values for the GraphViz layouter
 * plug-in.
 * 
 * @author ars
 */
public class GraphvizLayouterPreferenceInitializer extends AbstractPreferenceInitializer {

    @Override
    public void initializeDefaultPreferences() {
        IPreferenceStore store = GraphvizLayouterPreferenceStore.getInstance().getPreferenceStore();

        store.setDefault(GraphvizTool.PREF_TIMEOUT, GraphvizTool.PROCESS_DEF_TIMEOUT.get());
        store.setDefault(GraphvizLayoutProvider.PREF_GRAPHVIZ_REUSE_PROCESS,
                GraphvizLayoutProvider.REUSE_PROCESS_DEFAULT);
    }
}
