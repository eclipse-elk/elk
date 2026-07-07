/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.disco.debug;

import org.eclipse.elk.alg.disco.debug.views.DisCoDebugView;
import org.eclipse.elk.core.service.ILayoutListener;
import org.eclipse.elk.core.service.LayoutConnectorsService;
import org.eclipse.elk.core.service.LayoutMapping;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.ui.plugin.AbstractUIPlugin;
import org.osgi.framework.BundleContext;

/**
 * The activator class controls the plug-in life cycle.
 * 
 * This class started as a copy of DisCoDebugPlugin.java, commit 71bb8c2f542,
 * 2016-05-07. Changes are commented accordingly.
 */
public class DisCoDebugPlugin extends AbstractUIPlugin {

    /** The plug-in ID. */
    public static final String PLUGIN_ID = "org.eclipse.elk.alg.disco.debug"; //$NON-NLS-1$
    /** The shared instance. */
    private static DisCoDebugPlugin plugin;

    /**
     * The layout listener we will be using to update the view we are contributing.
     */
    private ILayoutListener layoutListener = new ILayoutListener() {
        @Override
        public void layoutAboutToStart(final LayoutMapping mapping, final IElkProgressMonitor progressMonitor) {
        }

        @Override
        public void layoutDone(final LayoutMapping mapping, final IElkProgressMonitor progressMonitor) {
            DisCoDebugView.updateWithGraph(mapping.getLayoutGraph());
        }
    };

    /**
     * The constructor.
     */
    public DisCoDebugPlugin() {
    }

    /*
     * (non-Javadoc)
     * 
     * @see org.eclipse.ui.plugin.AbstractUIPlugin#start(org.osgi.framework.BundleContext)
     */
    @Override
    public void start(final BundleContext context) throws Exception {
        super.start(context);
        plugin = this;

        LayoutConnectorsService.getInstance().addLayoutListener(layoutListener);
    }

    /*
     * (non-Javadoc)
     * 
     * @see org.eclipse.ui.plugin.AbstractUIPlugin#stop(org.osgi.framework.BundleContext)
     */
    @Override
    public void stop(final BundleContext context) throws Exception {
        plugin = null;
        super.stop(context);

        LayoutConnectorsService.getInstance().removeLayoutListener(layoutListener);
    }

    /**
     * Returns the shared instance.
     *
     * @return the shared instance
     */
    public static DisCoDebugPlugin getDefault() {
        return plugin;
    }

}
