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
package org.eclipse.elk.core.ui;

import org.eclipse.jface.resource.ImageDescriptor;
import org.eclipse.jface.resource.ImageRegistry;
import org.eclipse.ui.plugin.AbstractUIPlugin;
import org.osgi.framework.BundleContext;

/**
 * The activator class controls the plug-in life cycle.
 * 
 * @author haf
 * @author msp
 */
public class ElkUiPlugin extends AbstractUIPlugin {

    /** The plug-in ID. */
    public static final String PLUGIN_ID = "org.eclipse.elk.core.ui";
    
    /** Image ID for choice properties. */
    public static final String IMG_CHOICE = "obj16/prop_choice.gif";
    /** Image ID for double properties. */
    public static final String IMG_DOUBLE = "obj16/prop_double.gif";
    /** Image ID for integer properties. */
    public static final String IMG_INT = "obj16/prop_int.gif";
    /** Image ID for text properties. */
    public static final String IMG_TEXT = "obj16/prop_text.gif";
    /** Image ID for true Boolean properties. */
    public static final String IMG_TRUE = "obj16/prop_true.gif";
    /** Image ID for false Boolean properties. */
    public static final String IMG_FALSE = "obj16/prop_false.gif";

    /** the shared instance. */
    private static ElkUiPlugin plugin;

    /**
     * Returns the shared instance.
     * 
     * @return the shared instance
     */
    public static ElkUiPlugin getInstance() {
        return plugin;
    }

    @Override
    public void start(final BundleContext context) throws Exception {
        super.start(context);
        plugin = this;
    }
    
    @Override
    protected void initializeImageRegistry(final ImageRegistry reg) {
        for (String img : new String[] {IMG_CHOICE, IMG_DOUBLE, IMG_INT, IMG_TEXT, IMG_TRUE, IMG_FALSE}) {
            reg.put(img, getImageDescriptor("icons/" + img));
        }
    }

    @Override
    public void stop(final BundleContext context) throws Exception {
        plugin = null;
        super.stop(context);
    }

    /**
     * Returns an image descriptor for the image file at the given plug-in
     * relative path.
     * 
     * @param path the path
     * @return the image descriptor
     */
    public static ImageDescriptor getImageDescriptor(final String path) {
        return imageDescriptorFromPlugin(PLUGIN_ID, path);
    }
}
