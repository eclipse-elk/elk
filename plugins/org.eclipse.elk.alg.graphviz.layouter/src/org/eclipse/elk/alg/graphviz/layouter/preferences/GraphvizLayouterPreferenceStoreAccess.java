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

/**
 * Allows access to the {@link GraphvizLayouterPreferenceStore} in a ui save way. That is, the
 * eclipse.ui dependency is optional for this plugin. The methods of this class regard this by
 * returning default values if no preference store is available.
 * 
 * @author uru
 */
public final class GraphvizLayouterPreferenceStoreAccess {

    private static boolean uiAvailable = true;

    private GraphvizLayouterPreferenceStoreAccess() {
    }

    /**
     * @param name
     *            name of the requested preference
     * @return See {@link org.eclipse.jface.preference.IPreferenceStore#getInt(String)
     *         IPreferenceStore#getInt(String)}.
     */
    public static int getUISaveInt(final String name) {
        if (uiAvailable) {
            try {
                return GraphvizLayouterPreferenceStore.getInstance().getPreferenceStore()
                        .getInt(name);
            } catch (NoClassDefFoundError e) {
                uiAvailable = false;
            }
        }
        return 0;
    }

    /**
     * @param name
     *            name of the requested preference
     * @return See {@link org.eclipse.jface.preference.IPreferenceStore#getString(String)
     *         IPreferenceStore#getString(String)}.
     */
    public static String getUISaveString(final String name) {
        if (uiAvailable) {
            try {
                return GraphvizLayouterPreferenceStore.getInstance().getPreferenceStore()
                        .getString(name);
            } catch (NoClassDefFoundError e) {
                uiAvailable = false;
            }
        }
        return "";
    }

    /**
     * @param name
     *            name of the requested preference
     * @param ifFailed
     *            the boolean value to be returned if the preference store is not available.
     * @return See {@link org.eclipse.jface.preference.IPreferenceStore#getBoolean(String)
     *         IPreferenceStore#getBoolean(String)}.
     */
    public static boolean getUISaveBoolean(final String name, final boolean ifFailed) {
        if (uiAvailable) {
            try {
                return GraphvizLayouterPreferenceStore.getInstance().getPreferenceStore()
                        .getBoolean(name);
            } catch (NoClassDefFoundError e) {
                uiAvailable = false;
            }
        }
        return ifFailed;
    }
}
