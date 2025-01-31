/*******************************************************************************
 * Copyright (c) 2008, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.core.service;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.eclipse.core.runtime.Plugin;
import org.eclipse.core.runtime.preferences.InstanceScope;
import org.eclipse.elk.core.data.ILayoutMetaDataProvider;
import org.eclipse.elk.core.data.LayoutMetaDataService;
import org.eclipse.elk.core.service.util.MonitoredOperation;
import org.eclipse.elk.core.util.Pair;
import org.eclipse.jface.preference.IPreferenceStore;
import org.eclipse.ui.IWorkbenchPart;
import org.eclipse.ui.preferences.ScopedPreferenceStore;
import org.osgi.framework.Bundle;
import org.osgi.framework.BundleContext;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

/**
 * The activator class controls the plug-in life cycle.
 */
public final class ElkServicePlugin extends Plugin {

    /** The plug-in ID. */
    public static final String PLUGIN_ID = "org.eclipse.elk.core.service";

    /** The shared instance. */
    private static ElkServicePlugin plugin;
    
    /** The local preference store. */
    private IPreferenceStore preferenceStore;
    /** The executor service used to perform layout operations. */
    private ExecutorService executorService;
    /** map of currently running layout operations. */
    private final Multimap<Pair<IWorkbenchPart, Object>, MonitoredOperation> runningOperations = HashMultimap.create();
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Instance

    /**
     * Returns the shared instance.
     * 
     * @return the shared instance
     */
    public static ElkServicePlugin getInstance() {
        return plugin;
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Plug-in Lifecycle
    
    @Override
    public void start(final BundleContext context) throws Exception {
        super.start(context);
        plugin = this;
        
        // Load layout providers as a service
        loadLayoutProvidersThroughJavaServices();
    }

    @Override
    public void stop(final BundleContext context) throws Exception {
        LayoutMetaDataService.unload();
        LayoutConnectorsService.unload();
        if (executorService != null) {
            executorService.shutdown();
            executorService = null;
        }
        
        preferenceStore = null;
        plugin = null;
        super.stop(context);
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Accessors
    
    /**
     * Return the preference store associated with this plug-in.
     * 
     * Implementation is inspired by {@link org.eclipse.ui.plugin.AbstractUIPlugin#getPreferenceStore()}.
     */
    public IPreferenceStore getPreferenceStore() {
        // Create the preference store lazily.
        if (preferenceStore == null) {
            preferenceStore = new ScopedPreferenceStore(InstanceScope.INSTANCE, getBundle().getSymbolicName());

        }
        return preferenceStore;
    }
    
    /**
     * Return the executor service to use for automatic layout operations.
     */
    public synchronized ExecutorService getExecutorService() {
        if (executorService == null && plugin != null) {
            executorService = Executors.newCachedThreadPool();
        }
        return executorService;
    }
    
    /**
     * Return the map of currently running layout operations.
     */
    public Multimap<Pair<IWorkbenchPart, Object>, MonitoredOperation> getRunningOperations() {
        return runningOperations;
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Service and Extension Point Loading
    
    /**
     * Creates a new instance, loading the extension point information in the process.
     */
    @SuppressWarnings("unchecked")
    private void loadLayoutProvidersThroughJavaServices() {
        Set<Class<? extends ILayoutMetaDataProvider>> providerClasses = new HashSet<>();
        
        // Gather the classes of all ILayoutMetaDataProvider implementations
        for (Bundle bundle : getBundle().getBundleContext().getBundles()) {
            try {
                URL providerFile = bundle.getResource(
                        "/META-INF/services/" + ILayoutMetaDataProvider.class.getCanonicalName());
                if (providerFile != null) {
                    BufferedReader reader = new BufferedReader(new InputStreamReader(providerFile.openStream()));
                    
                    // Read the class names of all implementations that we can find
                    String line = reader.readLine();
                    while (line != null) {
                        providerClasses.add((Class<? extends ILayoutMetaDataProvider>) bundle.loadClass(line));
                        line = reader.readLine();
                    }
                }
            } catch (Exception e) {
                
            }
        }
        
        // Try instantiating and registering all of those classes
        LayoutMetaDataService service = LayoutMetaDataService.getInstance();
        for (Class<? extends ILayoutMetaDataProvider> providerClass : providerClasses) {
            ILayoutMetaDataProvider provider;
            try {
                provider = providerClass.getConstructor().newInstance();
                service.registerLayoutMetaDataProviders(provider);
            } catch (Exception e) {
            }
        }
    }
}
