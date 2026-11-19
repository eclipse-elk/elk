/*******************************************************************************
 * Copyright (c) 2016, 2019 Kiel University and others.
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
package org.eclipse.elk.core.debug;

import org.eclipse.elk.core.debug.model.ExecutionInfo;
import org.eclipse.elk.core.debug.model.ExecutionInfoModel;
import org.eclipse.elk.core.service.ILayoutListener;
import org.eclipse.elk.core.service.LayoutConnectorsService;
import org.eclipse.elk.core.service.LayoutMapping;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.ui.plugin.AbstractUIPlugin;
import org.osgi.framework.BundleContext;

/**
 * The activator class controls the plug-in life cycle
 */
public class ElkDebugPlugin extends AbstractUIPlugin {

	/** The plug-in ID. */
	public static final String PLUGIN_ID = "org.eclipse.elk.core.debug"; //$NON-NLS-1$
	/** The shared instance. */
	private static ElkDebugPlugin plugin;
	
	/** Our central execution info model that feeds our debug views. */
	private ExecutionInfoModel model = new ExecutionInfoModel();
	/** The layout listener we will be using to update the view we are contributing. */
	private ILayoutListener layoutListener = new ILayoutListener() {
        @Override
        public void layoutAboutToStart(final LayoutMapping mapping, final IElkProgressMonitor progressMonitor) {
        }
        
        @Override
        public void layoutDone(final LayoutMapping mapping, final IElkProgressMonitor progressMonitor) {
            model.addExecution(ExecutionInfo.fromProgressMonitor(progressMonitor));
        }
    };
	
    @Override
	public void start(BundleContext context) throws Exception {
		super.start(context);
		plugin = this;

        LayoutConnectorsService.getInstance().addLayoutListener(layoutListener);
	}
    
    @Override
	public void stop(BundleContext context) throws Exception {
		plugin = null;
		super.stop(context);
		
		LayoutConnectorsService.getInstance().removeLayoutListener(layoutListener);
	}

	/**
	 * Returns the shared instance
	 *
	 * @return the shared instance
	 */
	public static ElkDebugPlugin getDefault() {
		return plugin;
	}
	
	/**
	 * Returns our central execution info model.
	 */
	public ExecutionInfoModel getModel() {
	    return model;
	}

}
