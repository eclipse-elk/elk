/*******************************************************************************
 * Copyright (c) 2011, 2015 Kiel University and others.
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
package org.eclipse.elk.core.service;

import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;

import org.eclipse.core.runtime.CoreException;
import org.eclipse.core.runtime.IConfigurationElement;
import org.eclipse.core.runtime.Platform;
import org.eclipse.elk.core.service.internal.DefaultModule;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.ui.IWorkbenchPart;
import org.eclipse.ui.statushandlers.StatusManager;

import com.google.inject.Injector;

/**
 * A service class for layout connectors, which are provided by {@link ILayoutSetup} classes
 * registered through the extension point.
 *
 * @author msp
 */
public class LayoutConnectorsService {

    /** Identifier of the extension point for layout connectors. */
    protected static final String EXTP_ID_LAYOUT_CONNECTORS = "org.eclipse.elk.core.service.layoutConnectors";
    /** Name of the 'setup' element in the 'layoutConnectors' extension point. */
    protected static final String ELEMENT_SETUP = "setup";
    /** Name of the 'class' attribute in the extension points. */
    protected static final String ATTRIBUTE_CLASS = "class";
    /** Name of the 'priority' attribute in the extension points. */
    protected static final String ATTRIBUTE_PRIORITY = "priority";
    
    /** The singleton instance of the connectors service. */
    private static LayoutConnectorsService instance;
    
    /**
     * Returns the singleton instance of the connectors service.
     * 
     * @return the singleton instance
     */
    public static synchronized LayoutConnectorsService getInstance() {
        if (instance == null) {
            instance = new LayoutConnectorsService();
        }
        return instance;
    }
    
    /**
     * Unload any created instance in order to feed the garbage collector.
     */
    protected static synchronized void unload() {
        if (instance != null) {
            instance = null;
        }
    }
    
    /**
     * An entry read from the extension point.
     */
    private static class SetupEntry {
        /** The priority of this setup. */
        private int priority;
        /** The setup implementation. */
        private ILayoutSetup setup;
        /** The cached injector instance obtained from the setup. */
        private Injector injector;
    }
    
    
    /** List of registered layout setups. */
    private final List<SetupEntry> entries = new LinkedList<SetupEntry>();
    
    /** list of registered layout listeners. */
    private final List<ILayoutListener> layoutListeners = new LinkedList<ILayoutListener>();
    
    /**
     * Load all registered extensions for the layout connectors extension point.
     */
    public LayoutConnectorsService() {
        if (Platform.isRunning())
            loadLayoutSetupExtensions();
    }

    /**
     * Returns the most suitable layout setup injector for the given workbench and diagram part.
     * An implementation of {@link IDiagramLayoutConnector} can be obtained from such an injector.
     * 
     * @param workbenchPart the workbench part for which the injector should be fetched, or {@code null}
     * @param diagramPart the diagram part for which the injector should be fetched, or {@code null}
     * @return the most suitable injector, or {@code null} if none applies
     */
    public final Injector getInjector(final IWorkbenchPart workbenchPart, final Object diagramPart) {
        for (SetupEntry entry : entries) {
            if (workbenchPart == null) {
                if (entry.setup.supports(diagramPart)) {
                    return getInjector(entry);
                }
            } else if (entry.setup.supports(workbenchPart)
                    && (diagramPart == null || entry.setup.supports(diagramPart))) {
                return getInjector(entry);
            }
        }
        return null;
    }
    
    private Injector getInjector(final SetupEntry entry) {
        if (entry.injector == null) {
            entry.injector = entry.setup.createInjector(new DefaultModule());
        }
        return entry.injector;
    }
    
    /**
     * Returns the most suitable diagram layout connector for the given workbench and diagram part.
     * 
     * @param workbenchPart the workbench part for which the connector should be fetched, or {@code null}
     * @param diagramPart the diagram part for which the connector should be fetched, or {@code null}
     * @return the most suitable connector, or {@code null} if none applies
     */
    public final IDiagramLayoutConnector getConnector(final IWorkbenchPart workbenchPart,
            final Object diagramPart) {
        Injector injector = getInjector(workbenchPart, diagramPart);
        if (injector != null) {
            return injector.getInstance(IDiagramLayoutConnector.class);
        }
        return null;
    }
    
    /**
     * Add the given instance to the list of layout listeners.
     */
    public void addLayoutListener(final ILayoutListener listener) {
        layoutListeners.add(listener);
    }

    /**
     * Remove the given instance from the list of layout listeners.
     */
    public void removeLayoutListener(final ILayoutListener listener) {
        layoutListeners.remove(listener);
    }
    
    /**
     * Called by the {@link DiagramLayoutEngine} when automatic layout ist about to start.
     */
    protected void fireLayoutAboutToStart(final LayoutMapping mapping, final IElkProgressMonitor progressMonitor) {
        for (ILayoutListener listener : layoutListeners) {
            listener.layoutAboutToStart(mapping, progressMonitor);
        }
    }
    
    /**
     * Called by the {@link DiagramLayoutEngine} when automatic layout has been done.
     */
    protected void fireLayoutDone(final LayoutMapping mapping, final IElkProgressMonitor progressMonitor) {
        for (ILayoutListener listener : layoutListeners) {
            listener.layoutDone(mapping, progressMonitor);
        }
    }

    /**
     * Loads all layout setup extensions from the extension point.
     */
    private void loadLayoutSetupExtensions() {
        IConfigurationElement[] extensions = Platform.getExtensionRegistry()
                .getConfigurationElementsFor(EXTP_ID_LAYOUT_CONNECTORS);
        
        for (IConfigurationElement element : extensions) {
            try {
                if (ELEMENT_SETUP.equals(element.getName())) {
                    ILayoutSetup setup = (ILayoutSetup) element.createExecutableExtension(ATTRIBUTE_CLASS);
                    if (setup != null) {
                        SetupEntry entry = new SetupEntry();
                        entry.setup = setup;
                        String prioEntry = element.getAttribute(ATTRIBUTE_PRIORITY);
                        if (prioEntry != null) {
                            try {
                                entry.priority = Integer.parseInt(prioEntry);
                            } catch (NumberFormatException exception) {
                                // ignore exception
                            }
                        }
                        insertSorted(entry);
                    }
                    
                }
            } catch (CoreException exception) {
                StatusManager.getManager().handle(exception, ElkServicePlugin.PLUGIN_ID);
            }
        }
    }
    
    /**
     * Insert the given entry into the sorted list.
     */
    private void insertSorted(final SetupEntry entry) {
        ListIterator<SetupEntry> iter = entries.listIterator();
        while (iter.hasNext()) {
            SetupEntry next = iter.next();
            if (next.priority <= entry.priority) {
                iter.previous();
                break;
            }
        }
        iter.add(entry);
    }

}
