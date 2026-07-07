/*******************************************************************************
 * Copyright (c) 2011, 2019 Kiel University and others.
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

import java.io.IOException;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.concurrent.ExecutorService;

import org.eclipse.core.runtime.IStatus;
import org.eclipse.core.runtime.Status;
import org.eclipse.elk.core.IGraphLayoutEngine;
import org.eclipse.elk.core.LayoutConfigurator;
import org.eclipse.elk.core.data.LayoutAlgorithmResolver;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.SizeConstraint;
import org.eclipse.elk.core.service.util.MonitoredOperation;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.core.util.IElkCancelIndicator;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.core.util.IGraphElementVisitor;
import org.eclipse.elk.core.util.Maybe;
import org.eclipse.elk.core.util.Pair;
import org.eclipse.elk.core.validation.GraphValidator;
import org.eclipse.elk.core.validation.LayoutOptionValidator;
import org.eclipse.elk.graph.ElkGraphElement;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.properties.IPropertyHolder;
import org.eclipse.elk.graph.properties.MapPropertyHolder;
import org.eclipse.elk.graph.properties.Property;
import org.eclipse.emf.common.util.URI;
import org.eclipse.emf.ecore.resource.Resource;
import org.eclipse.emf.ecore.resource.ResourceSet;
import org.eclipse.emf.ecore.resource.impl.ResourceSetImpl;
import org.eclipse.jface.preference.IPreferenceStore;
import org.eclipse.ui.IWorkbenchPart;
import org.eclipse.ui.statushandlers.StatusManager;

import com.google.common.collect.Multimap;
import com.google.inject.ConfigurationException;
import com.google.inject.Inject;
import com.google.inject.Injector;
import com.google.inject.Provider;
import com.google.inject.Singleton;

/**
 * The entry class for automatic layout of graphical diagrams.
 * Use this class to perform automatic layout on the content of a workbench part that contains
 * a graph-based diagram. The mapping between the diagram and the layout graph structure is managed
 * by a {@link IDiagramLayoutConnector} implementation, which has to be registered using the
 * {@code layoutConnectors} extension point.
 * 
 * <p>Subclasses of this class can be bound in an {@link ILayoutSetup} injector for customization.
 * Note that this class is marked as {@link Singleton}, which means that exactly one instance is
 * created for each injector, i.e. for each registered {@link ILayoutSetup}.</p>
 */
@Singleton
public class DiagramLayoutEngine {
    
    /**
     * Configuration class for invoking the {@link DiagramLayoutEngine}.
     * Use a {@link LayoutConfigurator} to configure layout options:
     * <pre>
     * DiagramLayoutEngine.Parameters params = new DiagramLayoutEngine.Parameters();
     * params.addLayoutRun().configure(KNode.class)
     *         .setProperty(LayoutOptions.ALGORITHM, "org.eclipse.elk.algorithm.layered")
     *         .setProperty(LayoutOptions.SPACING, 30.0f)
     *         .setProperty(LayoutOptions.ANIMATE, true);
     * DiagramLayoutEngine.invokeLayout(workbenchPart, diagramPart, params);
     * </pre>
     * If multiple configurators are given, the layout is computed multiple times:
     * once for each configurator. This behavior can be used to apply different layout algorithms
     * one after another, e.g. first a node placer algorithm and then an edge router algorithm.
     * Example:
     * <pre>
     * DiagramLayoutEngine.Parameters params = new DiagramLayoutEngine.Parameters();
     * params.addLayoutRun().configure(KNode.class)
     *         .setProperty(LayoutOptions.ALGORITHM, "org.eclipse.elk.force");
     * params.addLayoutRun().setClearLayout(true).configure(KNode.class)
     *         .setProperty(LayoutOptions.ALGORITHM, "org.eclipse.elk.layered");
     * DiagramLayoutEngine.invokeLayout(workbenchPart, diagramPart, params);
     * </pre>
     * <b>Note:</b> By using the {@link LayoutConfigurator} approach as shown above, the Layout
     * view does not have any access to the configured values and hence will not work correctly.
     * In order to support the Layout view, use the {@link ILayoutConfigurationStore} interface instead.
     */
    public static final class Parameters {
        
        private List<IGraphElementVisitor> configurators = new LinkedList<IGraphElementVisitor>();
        private MapPropertyHolder globalSettings = new MapPropertyHolder();
        private boolean overrideDiagramConfig = true;
        
        /**
         * Set whether to override the configuration from the {@link ILayoutConfigurationStore}
         * provided by the diagram layout manager with the configuration provided by this setup.
         * The default is {@code true}.
         */
        public Parameters setOverrideDiagramConfig(final boolean override) {
            this.overrideDiagramConfig = override;
            return this;
        }
        
        /**
         * Returns the global settings, i.e. options that are not processed by layout algorithms,
         * but by layout managers or the layout engine itself.
         */
        public IPropertyHolder getGlobalSettings() {
            return globalSettings;
        }
        
        /**
         * Add a layout run with the given configurator. Each invocation of this method corresponds
         * to a separate layout execution on the whole input graph. This can be used to apply
         * multiple layout algorithms one after another, where each algorithm execution can reuse
         * results from previous executions.
         * 
         * @return the given configurator
         */
        public IGraphElementVisitor addLayoutRun(final IGraphElementVisitor configurator) {
            configurators.add(configurator);
            if (configurator instanceof LayoutConfigurator) {
                ((LayoutConfigurator) configurator).addFilter(LayoutConfigurator.OPTION_TARGET_FILTER);
            }
            return configurator;
        }
        
        /**
         * Convenience method for {@code addLayout(new LayoutConfigurator())}.
         */
        public LayoutConfigurator addLayoutRun() {
            return (LayoutConfigurator) addLayoutRun(new LayoutConfigurator());
        }
    }
    
    /** preference identifier for whether logging is enabled on progress monitors. */
    public static final String PREF_DEBUG_LOGGING = "elk.debug.logs";
    /** preference identifier for whether logged data are to be persisted. */
    public static final String PREF_DEBUG_STORE = "elk.debug.store";
    /** preference identifier for execution time measurement. */
    public static final String PREF_DEBUG_EXEC_TIME = "elk.debug.exectime";
        
    /**
     * Property for the diagram layout connector used for automatic layout. This property is
     * attached to layout mappings created by the {@code layout} methods.
     */
    public static final IProperty<IDiagramLayoutConnector> MAPPING_CONNECTOR
            = new Property<IDiagramLayoutConnector>("layoutEngine.diagramLayoutConnector");
    /**
     * Property for the status result of automatic layout. This property is attached to layout
     * mappings created by the {@code invokeLayout} methods.
     */
    public static final IProperty<IStatus> MAPPING_STATUS
            = new Property<IStatus>("layoutEngine.status");

    /**
     * Perform layout on the given workbench part with the given global options.
     * 
     * @param workbenchPart
     *            the workbench part for which layout is performed
     * @param diagramPart
     *            the parent diagram part for which layout is performed, or {@code null} if the whole
     *            diagram shall be arranged
     * @param animate
     *            if true, animation is activated (if supported by the diagram connector)
     * @param progressBar
     *            if true, a progress bar is displayed
     * @param layoutAncestors
     *            if true, layout is not only performed for the selected diagram part, but also for
     *            its ancestors
     * @param zoomToFit
     *            if true, automatic zoom-to-fit is activated (if supported by the diagram connector)
     * @return the layout mapping used in this operation
     */
    public static LayoutMapping invokeLayout(final IWorkbenchPart workbenchPart, final Object diagramPart,
            final boolean animate, final boolean progressBar, final boolean layoutAncestors,
            final boolean zoomToFit) {
        Parameters params = new Parameters();
        params.getGlobalSettings()
            .setProperty(CoreOptions.ANIMATE, animate)
            .setProperty(CoreOptions.PROGRESS_BAR, progressBar)
            .setProperty(CoreOptions.LAYOUT_ANCESTORS, layoutAncestors)
            .setProperty(CoreOptions.ZOOM_TO_FIT, zoomToFit);
        return invokeLayout(workbenchPart, diagramPart, params);
    }
    
    /**
     * Perform layout on the given workbench part with the given setup.
     * 
     * @param workbenchPart
     *            the workbench part for which layout is performed
     * @param diagramPart
     *            the parent diagram part for which layout is performed, or {@code null} if the whole
     *            diagram shall be layouted
     * @param params
     *            layout parameters, or {@code null} to use default values
     * @return the layout mapping used in this operation
     */
    public static LayoutMapping invokeLayout(final IWorkbenchPart workbenchPart, final Object diagramPart,
            final Parameters params) {
        return invokeLayout(workbenchPart, diagramPart, (IElkCancelIndicator) null, params);
    }

    /**
     * Perform layout on the given workbench part and diagram part. This static method creates an instance
     * of {@link DiagramLayoutEngine} using a {@link LayoutConnectorsService} injector and delegates the
     * operation to that instance.
     * 
     * <p>Depending on the {@code cancelIndicator} argument, different methods of the created engine
     * may be used: either the one taking an {@link IElkCancelIndicator} as argument, or the one taking
     * an {@link IElkProgressMonitor} as argument (the latter inherits from the former).</p>
     * 
     * <p>{@code workbenchPart} and {@code diagramPart} must not be {@code null} at the same time.</p>
     * 
     * @param workbenchPart
     *            the workbench part for which layout is performed, or {@code null}
     * @param diagramPart
     *            the parent diagram part for which layout is performed, or {@code null} if the whole
     *            diagram shall be layouted
     * @param cancelIndicator
     *            an {@link IElkCancelIndicator} to be evaluated repeatedly during the layout operation,
     *            or an {@link IElkProgressMonitor} to which progress is reported, or {@code null}
     * @param params
     *            layout parameters, or {@code null} to use default values
     * @return the layout mapping used in this operation, or {@code null} if the workbench part and diagram
     *            part cannot be identified by the {@link LayoutConnectorsService}
     */
    public static LayoutMapping invokeLayout(final IWorkbenchPart workbenchPart, final Object diagramPart,
            final IElkCancelIndicator cancelIndicator, final Parameters params) {
        Injector injector = LayoutConnectorsService.getInstance().getInjector(workbenchPart, diagramPart);

        if (injector != null) {
            try {
                DiagramLayoutEngine engine = injector.getInstance(DiagramLayoutEngine.class);
                if (cancelIndicator instanceof IElkProgressMonitor) {
                    return engine.layout(workbenchPart, diagramPart, (IElkProgressMonitor) cancelIndicator, params);
                } else {
                    return engine.layout(workbenchPart, diagramPart, cancelIndicator, params);
                }
            } catch (ConfigurationException exception) {
                IStatus status = new Status(IStatus.ERROR, ElkServicePlugin.PLUGIN_ID,
                        workbenchPart == null
                        ? "The Guice configuration for the given selection is inconsistent."
                        : "The Guice configuration for " + workbenchPart.getTitle() + " is inconsistent.",
                        exception);
                StatusManager.getManager().handle(status, StatusManager.SHOW);
            }
        } else {
            IStatus status = new Status(IStatus.ERROR, ElkServicePlugin.PLUGIN_ID,
                    workbenchPart == null
                    ? "No layout connector is available for the given selection."
                    : "No layout connector is available for " + workbenchPart.getTitle() + ".");
            StatusManager.getManager().handle(status, StatusManager.SHOW);
        }
        return null;
    }
    
    
    
    //--------------------- NON-STATIC PART (customizable via dependency injection) ---------------------//
    
    /**
     * The diagram layout connector used to import the layout graph and apply the resulting layout.
     */
    @Inject
    private IDiagramLayoutConnector connector;
    
    /**
     * The layout configuration manager for handling layout options.
     */
    @Inject
    private LayoutConfigurationManager configManager;
    
    /**
     * The graph layout engine that does the actual automatic layout computation.
     */
    @Inject
    private IGraphLayoutEngine graphLayoutEngine;
    
    /**
     * Resolves layout algorithms for nodes in the graph.
     */
    @Inject
    private LayoutAlgorithmResolver algorithmResolver;
    
    /**
     * Provider for validators of layout graphs.
     */
    @Inject
    private Provider<GraphValidator> graphValidatorProvider;
    
    /**
     * Provider for validators of layout options.
     */
    @Inject
    private Provider<LayoutOptionValidator> layoutOptionValidatorProvider;
    
    /**
     * Perform layout on the given workbench part and diagram part. The layout operation is wrapped in a
     * {@link MonitoredOperation} in order to ensure that building the layout graph and applying the layout
     * is both done in the UI thread. In case of a problem a message is shown to the user with the default
     * {@link StatusManager}.
     * 
     * <p>{@code workbenchPart} and {@code diagramPart} must not be {@code null} at the same time.</p>
     * 
     * @param workbenchPart
     *            the workbench part for which layout is performed, or {@code null}
     * @param diagramPart
     *            the parent diagram part for which layout is performed, or {@code null} if the whole
     *            diagram shall be layouted
     * @param cancelIndicator
     *            an {@link IElkCancelIndicator} to be evaluated repeatedly during the layout operation,
     *            or {@code null}
     * @param params
     *            layout parameters, or {@code null} to use default values
     * @return the layout mapping used in this operation
     */
    public LayoutMapping layout(final IWorkbenchPart workbenchPart, final Object diagramPart,
            final IElkCancelIndicator cancelIndicator, final Parameters params) {
        if (workbenchPart == null && diagramPart == null) {
            throw new NullPointerException();
        }
        
        final Parameters finalParams = params != null ? params : new Parameters();
        final Maybe<LayoutMapping> layoutMapping = Maybe.create();
        final Pair<IWorkbenchPart, Object> target = Pair.of(workbenchPart, diagramPart);
        final ExecutorService executorService = ElkServicePlugin.getInstance().getExecutorService();
        final MonitoredOperation monitoredOperation = new MonitoredOperation(executorService, cancelIndicator) {
            
            // First phase: build the layout graph
            @Override
            protected void preUIexec() {
                boolean layoutAncestors = finalParams.getGlobalSettings().getProperty(CoreOptions.LAYOUT_ANCESTORS);
                LayoutMapping mapping;
                if (layoutAncestors && workbenchPart != null) {
                    mapping = connector.buildLayoutGraph(workbenchPart, null);
                    mapping.setParentElement(diagramPart);
                } else {
                    mapping = connector.buildLayoutGraph(workbenchPart, diagramPart);
                }
                
                if (mapping != null && mapping.getLayoutGraph() != null) {
                    // Extract the diagram configuration
                    addDiagramConfig(finalParams, mapping);
                }
                
                layoutMapping.set(mapping);
            }

            // Second phase: execute layout algorithms
            @Override
            protected IStatus execute(final IElkProgressMonitor monitor) {
                if (monitor.isCanceled()) {
                    return Status.CANCEL_STATUS;
                }

                LayoutMapping mapping = layoutMapping.get();
                IStatus status;
                if (mapping != null && mapping.getLayoutGraph() != null) {
                    
                    // Perform the actual layout
                    status = layout(mapping, monitor, finalParams);
                    
                    // Stop earlier layout operations that are still running
                    if (!monitor.isCanceled()) {
                        stopEarlierOperations(target, getTimestamp());
                    }
                    
                    // Be sure that the monitor is closed
                    if (monitor.isRunning()) {
                        monitor.done();
                    }
                } else {
                    status = new Status(Status.WARNING, ElkServicePlugin.PLUGIN_ID,
                            "Unable to build the layout graph from the given selection.");
                }
                
                return status;
            }

            // Third phase: apply layout with animation
            @Override
            protected void postUIexec() {
                if (layoutMapping.get() != null) {
                    connector.applyLayout(layoutMapping.get(), finalParams.getGlobalSettings());
                }
            }
        };
        
        Multimap<Pair<IWorkbenchPart, Object>, MonitoredOperation> runningOperations = ElkServicePlugin
                .getInstance().getRunningOperations();
        synchronized (runningOperations) {
            runningOperations.put(target, monitoredOperation);
        }

        try {
            boolean progressBar = finalParams.getGlobalSettings().getProperty(CoreOptions.PROGRESS_BAR);
            if (progressBar) {
                // Perform layout with a progress bar
                monitoredOperation.runMonitored();
            } else {
                // Perform layout without a progress bar
                monitoredOperation.runUnmonitored();
            }
        } finally {
            synchronized (runningOperations) {
                runningOperations.remove(target, monitoredOperation);
            }
        }
        
        return layoutMapping.get();
    }
    
    /**
     * Stop all running operations whose timestamp is earlier than the given one.
     * 
     * @param target the layout target
     * @param time operations with a timestamp that is less than this are stopped
     */
    protected void stopEarlierOperations(final Pair<IWorkbenchPart, Object> target, final long time) {
        Multimap<Pair<IWorkbenchPart, Object>, MonitoredOperation> runningOperations = ElkServicePlugin
                .getInstance().getRunningOperations();
        synchronized (runningOperations) {
            for (MonitoredOperation operation : runningOperations.get(target)) {
                if (operation.getTimestamp() < time) {
                    operation.cancel();
                }
            }
        }
    }
    
    /**
     * Perform layout on the given workbench part and diagram part with a given progress monitor.
     * The three steps of the layout operation (build layout graph, invoke algorithms, apply layout)
     * are all executed in the same thread that calls this method. In case of a problem the resulting
     * {@link IStatus} is attached to the returned {@link LayoutMapping}, but is not reported to the user.
     * 
     * <p>{@code workbenchPart} and {@code diagramPart} must not be {@code null} at the same time.</p>
     * 
     * @param workbenchPart
     *            the workbench part for which layout is performed, or {@code null}
     * @param diagramPart
     *            the parent diagram part for which layout is performed, or {@code null} if the whole
     *            diagram shall be layouted
     * @param progressMonitor
     *            a progress monitor to which progress of the layout algorithm is reported,
     *            or {@code null} if no progress reporting is required
     * @param params
     *            layout parameters
     * @return the layout mapping used in this operation
     */
    public LayoutMapping layout(final IWorkbenchPart workbenchPart, final Object diagramPart,
            final IElkProgressMonitor progressMonitor, final Parameters params) {
        if (workbenchPart == null && diagramPart == null) {
            throw new NullPointerException();
        }
        
        final IElkProgressMonitor finalMonitor;
        if (progressMonitor == null) {
            IPreferenceStore prefStore = ElkServicePlugin.getInstance().getPreferenceStore();
            finalMonitor = new BasicProgressMonitor(0)
                    .withLogging(prefStore.getBoolean(PREF_DEBUG_LOGGING))
                    .withLogPersistence(prefStore.getBoolean(PREF_DEBUG_STORE))
                    .withExecutionTimeMeasurement(prefStore.getBoolean(PREF_DEBUG_EXEC_TIME));
        } else {
            finalMonitor = progressMonitor;
        }
        // SUPPRESS CHECKSTYLE NEXT MagicNumber
        finalMonitor.begin("Layout Diagram", 3);
        
        // Build the layout graph
        IElkProgressMonitor submon1 = finalMonitor.subTask(1);
        submon1.begin("Build layout graph", 1);
        LayoutMapping mapping = connector.buildLayoutGraph(workbenchPart, diagramPart);
        
        if (mapping != null && mapping.getLayoutGraph() != null) {
            // Extract the diagram configuration
            addDiagramConfig(params, mapping);
            submon1.done();
            
            // Perform the actual layout
            layout(mapping, finalMonitor.subTask(1), params);
            
            // Apply the layout to the diagram
            IElkProgressMonitor submon3 = finalMonitor.subTask(1);
            submon3.begin("Apply layout to the diagram", 1);
            connector.applyLayout(mapping, params.getGlobalSettings());
            submon3.done();
        } else {
            if (mapping == null) {
                mapping = new LayoutMapping(workbenchPart);
            }
            IStatus status = new Status(Status.WARNING, ElkServicePlugin.PLUGIN_ID,
                    "Unable to build the layout graph from the given selection.");
            mapping.setProperty(MAPPING_STATUS, status);
        }
        
        finalMonitor.done();
        return mapping;
    }
    
    /**
     * Create a diagram layout configuration and add it to the setup.
     */
    protected void addDiagramConfig(final Parameters params, final LayoutMapping layoutMapping) {
        LayoutConfigurator diagramConfig = configManager.createConfigurator(layoutMapping);
        if (params.configurators.isEmpty()) {
            params.addLayoutRun(diagramConfig);
        } else {
            ListIterator<IGraphElementVisitor> configIter = params.configurators.listIterator();
            while (configIter.hasNext()) {
                boolean isFirstConfig = !configIter.hasPrevious();
                IGraphElementVisitor setupConfig = configIter.next();
                if (setupConfig instanceof LayoutConfigurator) {
                    LayoutConfigurator layoutConfigurator = (LayoutConfigurator) setupConfig;
                    if (params.overrideDiagramConfig) {
                        if (isFirstConfig || layoutConfigurator.isClearLayout()) {
                            LayoutConfigurator newConfig;
                            if (configIter.hasNext()) {
                                newConfig = new LayoutConfigurator().overrideWith(diagramConfig);
                            } else {
                                newConfig = diagramConfig;
                            }
                            configIter.set(newConfig.overrideWith(layoutConfigurator));
                        }
                    } else {
                        layoutConfigurator.overrideWith(diagramConfig);
                    }
                }
            }
        }
    }
    
    /**
     * Handle the ancestors of the parent element if {@link CoreOptions#LAYOUT_ANCESTORS} is set.
     * For every ancestor node of the parent element (i.e. {@link LayoutMapping#getParentElement()}),
     * all containing elements that are not ancestors are excluded from layout.
     * 
     * @param mapping
     *            a mapping for the layout graph
     * @param params
     *            layout parameters
     */
    protected void handleAncestors(final LayoutMapping mapping, final Parameters params) {
        boolean layoutAncestors = params.getGlobalSettings().getProperty(CoreOptions.LAYOUT_ANCESTORS);
        if (layoutAncestors) {
            // Mark all parallel areas for exclusion from layout
            ElkGraphElement graphElem = mapping.getGraphMap().inverse().get(mapping.getParentElement());
            if (graphElem instanceof ElkNode && ((ElkNode) graphElem).getParent() != null) {
                if (params.configurators.isEmpty()) {
                    params.configurators.add(new LayoutConfigurator());
                }
                ElkNode node = (ElkNode) graphElem;
                do {
                    ElkNode parent = node.getParent();
                    for (ElkNode child : parent.getChildren()) {
                        if (child != node) {
                            for (IGraphElementVisitor c : params.configurators) {
                                if (c instanceof LayoutConfigurator) {
                                    IPropertyHolder childConfig = ((LayoutConfigurator) c).configure(child);
                                    // Do not layout the content of the child node
                                    childConfig.setProperty(CoreOptions.NO_LAYOUT, true);
                                    // Do not change the size of the child node
                                    childConfig.setProperty(CoreOptions.NODE_SIZE_CONSTRAINTS, SizeConstraint.fixed());
                                    // Do not move the ports of the child node
                                    childConfig.setProperty(CoreOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_POS);
                                }
                            }
                        }
                    }
                    node = parent;
                } while (node.getParent() != null);
            }
        }
    }
    
    /**
     * Perform layout on the given layout graph mapping. If zero or one layout configurator is
     * passed, the layout engine is executed exactly once. If multiple layout configurators are
     * passed, the layout engine is executed accordingly often, but the resulting layout is applied
     * only once. This is useful for composition of multiple algorithms that process only parts of
     * the graph. Layout listeners are notified before and after the layout has been computed.
     * 
     * @param mapping
     *            a mapping for the layout graph
     * @param progressMonitor
     *            a progress monitor to which progress of the layout algorithm is reported
     * @param params
     *            layout parameters
     * @return a status indicating success or failure
     */
    public IStatus layout(final LayoutMapping mapping, final IElkProgressMonitor progressMonitor,
            final Parameters params) {
        mapping.setProperty(MAPPING_CONNECTOR, connector);
        handleAncestors(mapping, params);
        
        LinkedList<IGraphElementVisitor> visitors = new LinkedList<IGraphElementVisitor>();
        visitors.add(algorithmResolver);
        
        // Set up graph validators
        if (params.getGlobalSettings().getProperty(CoreOptions.VALIDATE_OPTIONS)) {
            visitors.add(layoutOptionValidatorProvider.get());
        }
        if (params.getGlobalSettings().getProperty(CoreOptions.VALIDATE_GRAPH)) {
            visitors.add(graphValidatorProvider.get());
        }
        
        // Notify listeners of the to-be-executed layout
        LayoutConnectorsService.getInstance().fireLayoutAboutToStart(mapping, progressMonitor);
        
        IStatus status = null;
        if (params.configurators.isEmpty()) {
            // Perform layout without any extra configuration
            IGraphElementVisitor[] visitorsArray = visitors.toArray(new IGraphElementVisitor[visitors.size()]);
            status = layout(mapping, progressMonitor, visitorsArray);
        } else if (params.configurators.size() == 1) {
            // Perform layout once with an extra configuration
            visitors.addFirst(params.configurators.get(0));
            IGraphElementVisitor[] visitorsArray = visitors.toArray(new IGraphElementVisitor[visitors.size()]);
            status = layout(mapping, progressMonitor, visitorsArray);
        } else {
            // Perform layout multiple times with different configurations
            progressMonitor.begin("Diagram layout engine", params.configurators.size());
            ListIterator<IGraphElementVisitor> configIter = params.configurators.listIterator();
            while (configIter.hasNext()) {
                visitors.addFirst(configIter.next());
                IGraphElementVisitor[] visitorsArray = visitors.toArray(new IGraphElementVisitor[visitors.size()]);
                status = layout(mapping, progressMonitor, visitorsArray);
                if (!status.isOK()) {
                    break;
                }
                visitors.removeFirst();
                
                // If an additional layout configurator is attached to the graph, consider it in the future
                LayoutConfigurator addConfig =
                        mapping.getLayoutGraph().getProperty(LayoutConfigurator.ADD_LAYOUT_CONFIG);
                
                if (addConfig != null) {
                    ListIterator<IGraphElementVisitor> configIter2 = params.configurators.listIterator(
                            configIter.nextIndex());
                    while (configIter2.hasNext()) {
                        IGraphElementVisitor c = configIter2.next();
                        if (c instanceof LayoutConfigurator) {
                            ((LayoutConfigurator) c).overrideWith(addConfig);
                        }
                    }
                }
            }
            progressMonitor.done();
            
            // Log the final result to be displayed in our debug views
            progressMonitor.logGraph(mapping.getLayoutGraph(), "Result");
        }

        mapping.setProperty(MAPPING_STATUS, status);
        // Notify listeners of the executed layout
        LayoutConnectorsService.getInstance().fireLayoutDone(mapping, progressMonitor);
        return status;
    }
    
    /**
     * Perform layout on the given layout graph mapping. Layout listeners are <em>not</em> notified
     * in this method.
     * 
     * @param mapping
     *            a mapping for the layout graph
     * @param progressMonitor
     *            a progress monitor to which progress of the layout algorithm is reported; if the
     *            given monitor is not yet started, it is started in this method
     * @param visitors
     *            an optional array of graph element visitors to apply
     * @return a status indicating success or failure
     */
    public IStatus layout(final LayoutMapping mapping,
            final IElkProgressMonitor progressMonitor, final IGraphElementVisitor... visitors) {

        if (progressMonitor.isCanceled()) {
            return Status.CANCEL_STATUS;
        }
        boolean newTask = progressMonitor.begin("Diagram layout engine", 1);
        
        try {
            // Configure the layout graph by applying the given visitors
            if (visitors.length > 0) {
                ElkUtil.applyVisitorsWithValidation(mapping.getLayoutGraph(), visitors);
            }
            
            // Export the layout graph for debugging (this only does things if debug mode is enabled)
            progressMonitor.logGraph(mapping.getLayoutGraph(), "input");

            // Perform layout on the layout graph
            graphLayoutEngine.layout(mapping.getLayoutGraph(), progressMonitor.subTask(1));
            
            if (newTask) {
                progressMonitor.done();
                
                // Log the final result to be displayed in our debug views
                progressMonitor.logGraph(mapping.getLayoutGraph(), "result");
            }
            if (progressMonitor.isCanceled()) {
                return Status.CANCEL_STATUS;
            }
            
            // return a positive status
            return Status.OK_STATUS;
            
        } catch (Throwable exception) {
            return new Status(IStatus.ERROR, ElkServicePlugin.PLUGIN_ID,
                    "Failed to perform diagram layout.", exception);
        }
    }
    
    /**
     * Export the given layout graph in KGraph format.
     * 
     * @param graph the parent node of the layout graph
     */
    protected void exportLayoutGraph(final ElkNode graph) {
        URI exportUri = getExportURI(graph);
        if (exportUri != null) {
            // save the graph to a file
            ResourceSet resourceSet = new ResourceSetImpl();
            Resource resource = resourceSet.createResource(exportUri);
            resource.getContents().add(graph);
            try {
                resource.save(Collections.emptyMap());
            } catch (IOException e) {
                // ignore the exception and abort the layout graph exporting
            }
        }
    }
    
    /**
     * Return a file URI to use for exporting graphs.
     * 
     * @param graph the parent node of the layout graph
     */
    protected URI getExportURI(final ElkNode graph) {
        String path = ElkUtil.debugFolderPath("diagram_layout_engine")
                + Integer.toHexString(graph.hashCode()) + ".elkg";
        return URI.createFileURI(path);
    }

}
