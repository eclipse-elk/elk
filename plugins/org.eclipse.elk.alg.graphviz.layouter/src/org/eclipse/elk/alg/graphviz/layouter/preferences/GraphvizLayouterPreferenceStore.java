/*******************************************************************************
 * Copyright (c) 2015 Kiel University and others.
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

import org.eclipse.core.runtime.preferences.InstanceScope;
import org.eclipse.elk.alg.graphviz.layouter.GraphvizLayouterPlugin;
import org.eclipse.jface.preference.IPreferenceStore;
import org.eclipse.ui.preferences.ScopedPreferenceStore;

/**
 * In order to avoid extending {@link org.eclipse.ui.plugin.AbstractUIPlugin AbstractUIPlugin} by
 * the {@link GraphvizLayouterPreferenceStore} plugin, we maintain our own preference store here.
 * Usually the preference store is accessed using
 * {@link org.eclipse.ui.plugin.AbstractUIPlugin#getPreferenceStore()
 * AbstractUIPlugin#getPreferenceStore()}. The initialization of the store here is taken from the
 * just mentioned ui plugin.
 * 
 * @author uru
 */
public final class GraphvizLayouterPreferenceStore {

    private IPreferenceStore preferenceStore;

    private static GraphvizLayouterPreferenceStore instance;

    private GraphvizLayouterPreferenceStore() {
    }

    /**
     * @return an instance of {@link OgdfPreferenceStore}.
     */
    public static GraphvizLayouterPreferenceStore getInstance() {
        if (instance == null) {
            instance = new GraphvizLayouterPreferenceStore();
        }
        return instance;
    }

    /**
     * Functionally equal to {@link org.eclipse.ui.plugin.AbstractUIPlugin#getPreferenceStore()}.
     * 
     * @return the preference store
     */
    @SuppressWarnings("deprecation")
    public IPreferenceStore getPreferenceStore() {
        if (preferenceStore == null) {
            preferenceStore =
                    new ScopedPreferenceStore(new InstanceScope(), GraphvizLayouterPlugin
                            .getDefault().getBundle().getSymbolicName());
        }
        return preferenceStore;
    }
}
