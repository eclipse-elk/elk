/*******************************************************************************
 * Copyright (c) 2018, 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.test.framework.io;

import java.io.File;

import com.google.common.base.Strings;

/**
 * The path to a resource which is part of an Eclipse bundle. Such resources are always in the main ELK repository, the
 * path of which is expected to be stored in the system property {@code ELK_REPO}. The actual location is described as
 * a bundle ID and a relative path in that bundle.
 */
public class BundleResourcePath extends AbstractPropertyDependentResourcePath {
    
    /** The system property or environment variable that contains our base path. */
    public static final String PATH_PROPERTY = "ELK_REPO";

    /**
     * Creates a new instance for a resource in the given bundle.
     * 
     * @param bundleId
     *            the bundle's ID. This does not need to include any path relative to the ELK repository's root. The
     *            class automatically determines whether to look for it in the {@code plugins} or in the {@code test}
     *            folder based on the bundle ID.
     * @param filePath
     *            path to the resource, relative to the bundle's root folder.
     */
    public BundleResourcePath(final String bundleId, final String filePath) {
        if (Strings.isNullOrEmpty(bundleId)) {
            throw new IllegalArgumentException("The bundle ID cannot be empty.");
        }
        
        if (Strings.isNullOrEmpty(filePath)) {
            throw new IllegalArgumentException("The file path cannot be empty.");
        }
        
        initialize(bundleIdToPath(bundleId), filePath);
    }
    
    private String bundleIdToPath(final String bundleId) {
        // Find the ELK repository path
        String path = basePathForProperty(PATH_PROPERTY);

        // Bundles with an ID ending with "test" are located in the test directory. The others are in the plugins
        // directory. Other directories are not relevant to the test framework.
        if (bundleId.endsWith("test")) {
            path += "test" + File.separator;
        } else {
            path += "plugins" + File.separator;
        }
        
        // At this point the path specifies the point at which the directory with the bundle should be located. The
        // directory containing the bundle has the bundleId as its name.
        return path + bundleId;
    }

}
