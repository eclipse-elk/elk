/*******************************************************************************
 * Copyright (c) 2009, 2015 Kiel University and others.
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

import java.util.MissingResourceException;
import java.util.ResourceBundle;

/**
 * String externalization class for the ELK UI plugin.
 * 
 * @author msp
 */
public final class Messages {

    /** the bundle name. */
    private static final String BUNDLE_NAME = "org.eclipse.elk.core.ui.messages"; //$NON-NLS-1$
    /** the resource bundle instance. */
    private static final ResourceBundle RESOURCE_BUNDLE = ResourceBundle.getBundle(BUNDLE_NAME);

    /**
     * Hidden constructor.
     */
    private Messages() {
    }
    
    /**
     * Returns the string associated with the given key.
     * 
     * @param key key to look up in the {@code messages.properties} file
     * @return the associated string
     */
    public static String getString(final String key) {
        try {
            return RESOURCE_BUNDLE.getString(key);
        } catch (MissingResourceException e) {
            return '!' + key + '!';
        }
    }
}
