/*******************************************************************************
 * Copyright (c) 2018, 2019 Kiel University and others.
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
package org.eclipse.elk.alg.test.framework.io;

/**
 * Subclasses of this class represent resource paths whose base path depends on a system property or environment
 * variable. Subclasses work just like direct subclasses of {@link AbstractResourcePath}, but can use the
 * {@link #basePathForProperty(String)} method to resolve the base path. 
 */
public abstract class AbstractPropertyDependentResourcePath extends AbstractResourcePath {
    
    /**
     * Determines the base path by querying the system property or environment variable with the given name. The
     * returned path will always end with a {@code /}.
     * 
     * @param propertyName
     *            the property name.
     * @return the base path.
     * @throws IllegalStateException
     *             if no system property and environment variable with the given name exists.
     */
    protected static String basePathForProperty(final String propertyName) {
        String path = System.getProperty(propertyName);
        if (path == null) {
            path = System.getenv(propertyName);
        }
        
        if (path == null) {
            throw new IllegalStateException("The system property or environment variable '"
                    + propertyName + "' needs to be set.");
        }
        
        return path;
    }
    
}
