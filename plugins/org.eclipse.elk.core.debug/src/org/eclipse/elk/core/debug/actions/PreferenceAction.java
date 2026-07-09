/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
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
package org.eclipse.elk.core.debug.actions;

import org.eclipse.elk.core.debug.ElkDebugPlugin;
import org.eclipse.elk.core.service.ElkServicePlugin;
import org.eclipse.jface.action.Action;
import org.eclipse.jface.preference.IPreferenceStore;
import org.eclipse.jface.util.IPropertyChangeListener;
import org.eclipse.jface.util.PropertyChangeEvent;

/**
 * Action that updates logging preferences. Call {@link #dispose()} after use.
 */
public class PreferenceAction extends Action {

    /** Identifier string for this action. */
    private static final String ACTION_ID = "org.eclipse.elk.debug.actions.preference";
    /** Relative path to the icon to use for this action. */
    private static final String ICON_PATH = "icons/preference.gif";
    
    /** Key of the preference we control. */
    private final String preferenceId;
    /** Listener for property changes to update our check state. */
    private final PropertyChangeListener propertyChangeListener = new PropertyChangeListener();
    
    /**
     * Creates an action that controls the preference with the given ID.
     */
    public PreferenceAction(final String preferenceId) {
        super("&Log", Action.AS_CHECK_BOX);
        
        setId(ACTION_ID);
        setToolTipText("Enables or disables logging information displayed in this view.");
        setImageDescriptor(ElkDebugPlugin.imageDescriptorFromPlugin(ElkDebugPlugin.PLUGIN_ID, ICON_PATH));
        
        this.preferenceId = preferenceId;
        updateCheckedState();
        
        // Listen for preference changes
        ElkServicePlugin.getInstance().getPreferenceStore().addPropertyChangeListener(propertyChangeListener);
    }
    
    /**
     * Shuts down this action.
     */
    public void dispose() {
        ElkServicePlugin.getInstance().getPreferenceStore().removePropertyChangeListener(propertyChangeListener);
    }
    
    /**
     * Updates the checked property of this action.
     */
    public void updateCheckedState() {
        setChecked(ElkServicePlugin.getInstance().getPreferenceStore().getBoolean(preferenceId));
    }

    @Override
    public void run() {
        IPreferenceStore prefStore = ElkServicePlugin.getInstance().getPreferenceStore();
        prefStore.setValue(preferenceId, !prefStore.getBoolean(preferenceId));
    }
    
    /**
     * Listens to changes to our property.
     */
    private class PropertyChangeListener implements IPropertyChangeListener {
        @Override
        public void propertyChange(PropertyChangeEvent event) {
            if (event.getProperty().equals(preferenceId)) {
                updateCheckedState();
            }
        }
    }

}
