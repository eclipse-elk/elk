/*******************************************************************************
 * Copyright (c) 2013, 2022 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.libavoid;

import org.eclipse.core.runtime.Plugin;
import org.eclipse.elk.alg.libavoid.server.LibavoidServerPool;
import org.osgi.framework.BundleContext;

/**
 * The activator class controls the plug-in life cycle.
 * 
 * @author uru
 */
public class LibavoidPlugin extends Plugin {

    /** the plug-in ID. */
    public static final String PLUGIN_ID = "org.adaptagrams.cola.libavoid"; //$NON-NLS-1$

    /** the shared instance. */
    private static LibavoidPlugin plugin;
    
    /**
     * The constructor.
     */
    public LibavoidPlugin() {
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void start(final BundleContext context) throws Exception {
        super.start(context);
        plugin = this;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void stop(final BundleContext context) throws Exception {
        plugin = null;
        LibavoidServerPool.INSTANCE.dispose();
        super.stop(context);
    }

    /**
     * Returns the shared instance.
     *
     * @return the shared instance
     */
    public static LibavoidPlugin getDefault() {
        return plugin;
    }

}
