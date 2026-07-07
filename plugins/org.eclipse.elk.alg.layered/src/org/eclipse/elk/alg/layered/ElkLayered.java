/*******************************************************************************
 * Copyright (c) 2010, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Deque;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.util.Set;

import org.eclipse.elk.alg.layered.components.ComponentsProcessor;
import org.eclipse.elk.alg.layered.compound.CompoundGraphPostprocessor;
import org.eclipse.elk.alg.layered.compound.CompoundGraphPreprocessor;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPadding;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.CrossingMinimizationStrategy;
import org.eclipse.elk.alg.layered.options.GraphProperties;
import org.eclipse.elk.alg.layered.options.GreedySwitchType;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.UnsupportedGraphException;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.ContentAlignment;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.options.SizeConstraint;
import org.eclipse.elk.core.options.SizeOptions;
import org.eclipse.elk.core.testing.TestController;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.core.util.Pair;

/**
 * The main entry point into ELK Layered. ELK Layered is a layout algorithm after the layered
 * layout method proposed by Sugiyama et al. It is structured into five main phases: cycle breaking,
 * layering, crossing minimization, node placement, and edge routing. Before these phases and after
 * the last phase so called intermediate layout processors can be inserted that do some kind of pre
 * or post processing. Implementations of the different main phases specify the intermediate layout
 * processors they require, which are automatically collected and inserted between the main phases.
 * The layout provider itself also specifies some dependencies.
 *
 * <pre>
 *           Intermediate Layout Processors
 * ---------------------------------------------------
 * |         |         |         |         |         |
 * |   ---   |   ---   |   ---   |   ---   |   ---   |
 * |   | |   |   | |   |   | |   |   | |   |   | |   |
 * |   | |   |   | |   |   | |   |   | |   |   | |   |
 *     | |       | |       | |       | |       | |
 *     | |       | |       | |       | |       | |
 *     ---       ---       ---       ---       ---
 *   Phase 1   Phase 2   Phase 3   Phase 4   Phase 5
 * </pre>
 *
 * <p>To use ELK Layered to layout a given graph, there are three possibilities depending on the kind
 * of graph that is to be laid out:</p>
 * <ol>
 *   <li>{@link #doLayout(LGraph, IElkProgressMonitor)} computes a layout for the given graph, without
 *     any subgraphs it might have.</li>
 *   <li>{@link #doCompoundLayout(LGraph, IElkProgressMonitor)} computes a layout for the given graph
 *     and for its subgraphs, if any. (Subgraphs are attached to nodes through the
 *     {@link InternalProperties#NESTED_LGRAPH} property.)</li>
 *   <li>If you have an {@code ElkNode} instead of an {@code LGraph}, you might want to use
 *     {@link LayeredLayoutProvider#layout(org.eclipse.elk.graph.ElkNode, IElkProgressMonitor)}
 *     instead.</li>
 * </ol>
 * <p>In addition to regular layout runs, this class provides methods for automatic unit testing based
 * around the concept of <em>test runs</em>. A test run is executed as follows:</p>
 * <ol>
 *   <li>Call {@link #prepareLayoutTest(LGraph)} to start a new run. The given graph might be split
 *       into its connected components, which are put into the returned state object.</li>
 *   <li>Call one of the actual test methods. {@link #runLayoutTestStep(TestExecutionState)} runs the
 *       next step of the algorithm. {@link #runLayoutTestUntil(Class, TestExecutionState)} runs the
 *       algorithm until a given layout processor has finished executing (its sibling,
 *       {@link #runLayoutTestUntil(Class, boolean, TestExecutionState)}, can also stop just before a
 *       given layout processor starts executing). All of these methods resume execution from where the
 *       algorithm has stopped previously.</li>
 * </ol>
 * <p>ELK Layered also supports the ELK unit test framework through being a white box testable algorithm.</p>
 *
 * @see ILayoutProcessor
 * @see GraphConfigurator
 * @see IHierarchyAwareLayoutProcessor
 *
 * @author msp
 * @author cds
 */
public final class ElkLayered {

    ////////////////////////////////////////////////////////////////////////////////
    // Variables

    /** the algorithm's current configuration. */
    private final GraphConfigurator graphConfigurator = new GraphConfigurator();
    /** connected components processor. */
    private final ComponentsProcessor componentsProcessor = new ComponentsProcessor();
    /** compound graph preprocessor. */
    private final CompoundGraphPreprocessor compoundGraphPreprocessor = new CompoundGraphPreprocessor();
    /** compound graph postprocessor. */
    private final CompoundGraphPostprocessor compoundGraphPostprocessor = new CompoundGraphPostprocessor();
    /** Test controller for a white box test. */
    private TestController testController = null;


    ////////////////////////////////////////////////////////////////////////////////
    // Regular Layout

    /**
     * Does a layout on the given graph. If the graph contains compound nodes (see class documentation),
     * the nested graphs are ignored.
     *
     * @param lgraph the graph to layout
     * @param monitor a progress monitor to show progress information in, or {@code null}
     * @see #doCompoundLayout(LGraph, IElkProgressMonitor)
     */
    public void doLayout(final LGraph lgraph, final IElkProgressMonitor monitor) {
        IElkProgressMonitor theMonitor = monitor;
        if (theMonitor == null) {
            theMonitor = new BasicProgressMonitor().withMaxHierarchyLevels(0);
        }
        theMonitor.begin("Layered layout", 1);

        // Update the modules depending on user options
        graphConfigurator.prepareGraphForLayout(lgraph);

        // Split the input graph into components and perform layout on them
        List<LGraph> components = componentsProcessor.split(lgraph);
        if (components.size() == 1) {
            // Execute layout on the sole component using the top-level progress monitor
            layout(components.get(0), theMonitor);
        } else {
            // Execute layout on each component using a progress monitor subtask
            float compWork = 1.0f / components.size();
            for (LGraph comp : components) {
                if (monitor.isCanceled()) {
                    return;
                }
                layout(comp, theMonitor.subTask(compWork));
            }
        }
        componentsProcessor.combine(components, lgraph);

        // Resize the resulting graph, according to minimal size constraints and such
        resizeGraph(lgraph);

        theMonitor.done();
    }


    ////////////////////////////////////////////////////////////////////////////////
    // Compound Graph Layout

    /**
     * Does a layout on the given compound graph. Connected components processing is currently not
     * supported.
     *
     * @param lgraph the graph to layout
     * @param monitor a progress monitor to show progress information in, or {@code null}
     */
    public void doCompoundLayout(final LGraph lgraph, final IElkProgressMonitor monitor) {
        IElkProgressMonitor theMonitor = monitor;
        if (theMonitor == null) {
            theMonitor = new BasicProgressMonitor().withMaxHierarchyLevels(0);
        }
        theMonitor.begin("Layered layout", 2); // SUPPRESS CHECKSTYLE MagicNumber

        // Preprocess the compound graph by splitting cross-hierarchy edges
        notifyProcessorReady(lgraph, compoundGraphPreprocessor);
        compoundGraphPreprocessor.process(lgraph, theMonitor.subTask(1));
        notifyProcessorFinished(lgraph, compoundGraphPreprocessor);

        hierarchicalLayout(lgraph, theMonitor.subTask(1));

        // Postprocess the compound graph by combining split cross-hierarchy edges
        notifyProcessorReady(lgraph, compoundGraphPostprocessor);
        compoundGraphPostprocessor.process(lgraph, theMonitor.subTask(1));
        notifyProcessorFinished(lgraph, compoundGraphPostprocessor);

        theMonitor.done();
    }

    /**
     * Processors can be marked as operating on the full hierarchy by using the {@link IHierarchyAwareLayoutProcessor}
     * interface.
     *
     * All graphs are collected using a breadth-first search and this list is reversed, so that for each graph, all
     * following graphs are on the same hierarchy level or higher, i.e. closer to the parent graph. Each graph then has
     * a unique configuration of ELK Layered, which is comprised of a sequence of processors. The processors can vary
     * depending on the characteristics of each graph. The list of graphs and their algorithms is then traversed. If a
     * processor is not hierarchical it is simply executed. If it it is hierarchical and this graph is not the root
     * graph, this processor is skipped and the algorithm is paused until the processor has been executed on the root
     * graph. Then the algorithm is continued, starting with the level lowest in the hierarchy, i.e. furthest away from
     * the root graph.
     */
    private void hierarchicalLayout(final LGraph lgraph, final IElkProgressMonitor monitor) {
        // Perform a reversed breadth first search: The graphs in the lowest hierarchy come first.
        Collection<LGraph> graphs = collectAllGraphsBottomUp(lgraph);

        // We have to make sure that hierarchical processors don't break the control flow 
        //  of the following layout execution (see e.g. #228). To be precise, if the root node of 
        //  the hierarchical graph doesn't include a hierarchical processor, nor may any of the children.
        reviewAndCorrectHierarchicalProcessors(lgraph, graphs);
        
        // Get list of processors for each graph, since they can be different.
        // Iterators are used, so that processing of a graph can be paused and continued easily.
        int work = 0;
        List<Pair<LGraph, Iterator<ILayoutProcessor<LGraph>>>> graphsAndAlgorithms = new ArrayList<>();
        for (LGraph g : graphs) {
            graphConfigurator.prepareGraphForLayout(g);
            List<ILayoutProcessor<LGraph>> processors = g.getProperty(InternalProperties.PROCESSORS);
            work += processors.size();
            Iterator<ILayoutProcessor<LGraph>> algorithm = processors.iterator();
            graphsAndAlgorithms.add(Pair.of(g, algorithm));
        }

        monitor.begin("Recursive hierarchical layout", work);

        // When the root graph has finished layout, the layout is complete.
        int slotIndex = 0;
        Iterator<ILayoutProcessor<LGraph>> rootProcessors = getProcessorsForRootGraph(graphsAndAlgorithms);
        while (rootProcessors.hasNext()) {
            // Layout from bottom up
            for (Pair<LGraph, Iterator<ILayoutProcessor<LGraph>>> graphAndAlgorithm : graphsAndAlgorithms) {
                Iterator<ILayoutProcessor<LGraph>> processors = graphAndAlgorithm.getSecond();
                LGraph graph = graphAndAlgorithm.getFirst();

                while (processors.hasNext()) {
                    ILayoutProcessor<LGraph> processor = processors.next();
                    if (!(processor instanceof IHierarchyAwareLayoutProcessor)) {
                        // Output debug graph
                        // elkjs-exclude-start
                        if (monitor.isLoggingEnabled()) {
                            DebugUtil.logDebugGraph(monitor, graph, slotIndex,
                                    "Before " + processor.getClass().getSimpleName());
                        }
                        // elkjs-exclude-end
                        
                        notifyProcessorReady(graph, processor);
                        processor.process(graph, monitor.subTask(1));
                        notifyProcessorFinished(graph, processor);
                        
                        slotIndex++;
                        
                    } else if (isRoot(graph)) {
                        // Output debug graph
                        // elkjs-exclude-start
                        if (monitor.isLoggingEnabled()) {
                            DebugUtil.logDebugGraph(monitor, graph, slotIndex,
                                    "Before " + processor.getClass().getSimpleName());
                        }
                        // elkjs-exclude-end
                        
                        // If processor operates on the full hierarchy, it must be executed on the root
                        notifyProcessorReady(graph, processor);
                        processor.process(graph, monitor.subTask(1));
                        notifyProcessorFinished(graph, processor);
                        
                        slotIndex++;
                        
                        // Continue operation with the graph at the bottom of the hierarchy
                        break;
                        
                    } else { // operates on full hierarchy and is not root graph
                        // skip this processor and pause execution until root graph has processed.
                        break;
                    }
                }
            }
        }

        // Graph debug output
        // elkjs-exclude-start
        if (monitor.isLoggingEnabled()) {
            DebugUtil.logDebugGraph(monitor, lgraph, slotIndex, "Finished");
        }
        // elkjs-exclude-end

        monitor.done();
    }

    /**
     * Implements a breadth first search in compound graphs with reversed order. This way the
     * innermost graphs are first in the list, followed by one level further up, etc.
     *
     * @param root
     *            the root graph
     * @return Graphs in breadth first search in compound graphs in reverse order.
     */
    private Collection<LGraph> collectAllGraphsBottomUp(final LGraph root) {
        Deque<LGraph> collectedGraphs = new ArrayDeque<>();
        Deque<LGraph> continueSearchingTheseGraphs = new ArrayDeque<>();
        collectedGraphs.push(root);
        continueSearchingTheseGraphs.push(root);

        while (!continueSearchingTheseGraphs.isEmpty()) {
            LGraph nextGraph = continueSearchingTheseGraphs.pop();
            for (LNode node : nextGraph.getLayerlessNodes()) {
                if (hasNestedGraph(node)) {
                    LGraph nestedGraph = node.getNestedGraph();
                    collectedGraphs.push(nestedGraph);
                    continueSearchingTheseGraphs.push(nestedGraph);
                }
            }
        }
        return collectedGraphs;
    }
    
    /**
     * It is not permitted that any of the child-graphs specifies a hierarchical 
     * layout processor ({@link IHierarchyAwareLayoutProcessor}) that is not specified by the root node.
     * 
     * It depends on the concrete processor how this is fixed.   
     * 
     * @param root the root graph
     * @param graphs all graphs of the handled hierarchy
     */
    private void reviewAndCorrectHierarchicalProcessors(final LGraph root, final Collection<LGraph> graphs) {
        // Crossing minimization
        //  overwrite invalid child configuration (only layer sweep is hierarchical)
        CrossingMinimizationStrategy parentCms = root.getProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY);
        graphs.forEach(child -> {
            CrossingMinimizationStrategy childCms = 
                    child.getProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY);
            if (childCms != parentCms) {
                throw new UnsupportedGraphException("The hierarchy aware processor " + childCms + " in child node "
                        + child + " is only allowed if the root node specifies the same hierarchical processor.");
            }
        });
        
        // Greedy switch (simply copy the behavior of the root to all children)
        final GreedySwitchType rootType =
                root.getProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_HIERARCHICAL_TYPE);
        graphs.forEach(
                g -> g.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_HIERARCHICAL_TYPE, rootType));
    }

    private Iterator<ILayoutProcessor<LGraph>> getProcessorsForRootGraph(
            final List<Pair<LGraph, Iterator<ILayoutProcessor<LGraph>>>> graphsAndAlgorithms) {
        return graphsAndAlgorithms.get(graphsAndAlgorithms.size() - 1).getSecond();
    }

    private boolean isRoot(final LGraph graph) {
        return graph.getParentNode() == null;
    }

    private boolean hasNestedGraph(final LNode node) {
        return node.getNestedGraph() != null;
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Layout Testing

    /**
     * The state of a test execution is held in an instance of this class.
     */
    public static class TestExecutionState {
        /** list of graphs that are currently being laid out. */
        private List<LGraph> graphs;
        /** index of the processor that is to be executed next during a layout test. */
        private int step;

        /**
         * Return the list of graphs that are currently being laid out.
         *
         * @return the graphs under test
         */
        public List<LGraph> getGraphs() {
            return graphs;
        }

        /**
         * Return the index of the processor that is to be executed next during a layout test.
         *
         * @return the index of the next step
         */
        public int getStep() {
            return step;
        }
    }

    /**
     * Prepares a test run of the layout algorithm. After this method has run, call
     * {@link #layoutTestStep()} as often as there are layout processors.
     *
     * @param lgraph the input graph to initialize the test run with.
     * @return the test execution state
     */
    public TestExecutionState prepareLayoutTest(final LGraph lgraph) {
        TestExecutionState state = new TestExecutionState();

        // update the modules depending on user options
        graphConfigurator.prepareGraphForLayout(lgraph);

        // split the input graph into components
        state.graphs = componentsProcessor.split(lgraph);

        return state;
    }

    /**
     * Checks if the current test run still has processors to be executed for the algorithm to finish.
     *
     * @param state the current test execution state
     * @return {@code true} if the current test run has not finished yet. If there is no current
     *         test run, the result is undefined.
     */
    public boolean isLayoutTestFinished(final TestExecutionState state) {
        LGraph graph = state.graphs.get(0);
        List<ILayoutProcessor<LGraph>> algorithm = graph.getProperty(InternalProperties.PROCESSORS);
        return algorithm != null && state.step >= algorithm.size();
    }

    /**
     * Runs the algorithm on the current test graphs up to the point where the given phase or
     * processor has finished executing. If parts of the algorithm were already executed using this
     * or other layout test methods, execution is resumed from there. If the given phase or
     * processor is not among those processors that have not yet executed, an exception is thrown.
     * Also, if there is no current layout test run, an exception is thrown.
     *
     * @param phase the phase or processor to stop after
     * @param inclusive {@code true} if the specified phase should be executed as well
     * @param state the current test execution state
     * @throws IllegalArgumentException
     *             if the given layout processor is not part of the processors that are still to be
     *             executed.
     */
    public void runLayoutTestUntil(final Class<? extends ILayoutProcessor<LGraph>> phase,
            final boolean inclusive, final TestExecutionState state) {

        List<ILayoutProcessor<LGraph>> algorithm = state.graphs.get(0).getProperty(InternalProperties.PROCESSORS);

        // check if the given phase exists in our current algorithm configuration
        boolean phaseExists = false;
        ListIterator<ILayoutProcessor<LGraph>> algorithmIterator = algorithm.listIterator(state.step);
        int phaseIndex = state.step;

        while (algorithmIterator.hasNext() && !phaseExists) {
            if (algorithmIterator.next().getClass().equals(phase)) {
                phaseExists = true;

                if (inclusive) {
                    phaseIndex++;
                }
            } else {
                phaseIndex++;
            }
        }

        if (!phaseExists) {
            // FIXME actually, we want to know when a processor is not
            //  part of the algorithm's configuration because this might be
            //  wrong behavior.
            // However, in the current test framework there is no way
            //  to differentiate between 'it's ok' and 'it's not'.
            // throw new IllegalArgumentException(
            // "Given processor not part of the remaining algorithm.");
            System.err.println("Given processor " + phase + " not part of the remaining algorithm.");
        }

        // perform the layout up to and including that phase
        algorithmIterator = algorithm.listIterator(state.step);
        for (; state.step < phaseIndex; state.step++) {
            layoutTest(state.graphs, algorithmIterator.next());
        }
    }

    /**
     * Performs the {@link #runLayoutTestUntil(Class, boolean)} methods with {@code inclusive} set
     * to {@code true}.
     *
     * @param phase the phase or processor to stop after
     * @param state the current test execution state
     * @see ElkLayered#runLayoutTestUntil(Class, boolean)
     */
    public void runLayoutTestUntil(final Class<? extends ILayoutProcessor<LGraph>> phase,
            final TestExecutionState state) {

        runLayoutTestUntil(phase, true, state);
    }

    /**
     * Runs the next step of the current layout test run. Throws exceptions if no layout test run is
     * currently active or if the current run has finished.
     *
     * @param state the current test execution state
     * @throws IllegalStateException if the given state has finished executing
     */
    public void runLayoutTestStep(final TestExecutionState state) {
        if (isLayoutTestFinished(state)) {
            throw new IllegalStateException("Current layout test run has finished.");
        }

        // perform the next layout step
        List<ILayoutProcessor<LGraph>> algorithm = state.graphs.get(0).getProperty(InternalProperties.PROCESSORS);
        layoutTest(state.graphs, algorithm.get(state.step));
        state.step++;
    }

    /**
     * Returns the current list of layout processors that make up the algorithm. This list is only
     * valid and meaningful while a layout test is being run.
     *
     * @param state the current test execution state
     * @return the algorithm's current configuration.
     */
    public List<ILayoutProcessor<LGraph>> getLayoutTestConfiguration(final TestExecutionState state) {
        return state.graphs.get(0).getProperty(InternalProperties.PROCESSORS);
    }
    
    /**
     * Installs the given test controller to be notified of layout events.
     * 
     * @param testController the test controller to be installed. May be {@code null}.
     */
    public void setTestController(final TestController testController) {
        this.testController = testController;
    }
    
    /**
     * Notifies the test controller (if installed) that the given processor is ready to start processing the given
     * graph. If the graph is the root graph, the corresponding notification is triggered as well.
     */
    private void notifyProcessorReady(final LGraph lgraph, final ILayoutProcessor<?> processor) {
        if (testController != null) {
            if (isRoot(lgraph)) {
                testController.notifyRootProcessorReady(lgraph, processor);
            } else {
                testController.notifyProcessorReady(lgraph, processor);
            }
        }
    }

    /**
     * Notifies the test controller (if installed) that the given processor has finished processing the given
     * graph. If the graph is the root graph, the corresponding notification is triggered as well.
     */
    private void notifyProcessorFinished(final LGraph lgraph, final ILayoutProcessor<?> processor) {
        if (testController != null) {
            if (isRoot(lgraph)) {
                testController.notifyRootProcessorFinished(lgraph, processor);
            } else {
                testController.notifyProcessorFinished(lgraph, processor);
            }
        }
    }


    ////////////////////////////////////////////////////////////////////////////////
    // Actual Layout

    /**
     * Perform the five phases of the layered layouter.
     *
     * @param lgraph the graph that is to be laid out
     * @param monitor a progress monitor
     */
    private void layout(final LGraph lgraph, final IElkProgressMonitor monitor) {
        boolean monitorWasAlreadyRunning = monitor.isRunning();
        if (!monitorWasAlreadyRunning) {
            monitor.begin("Component Layout", 1);
        }
        
        List<ILayoutProcessor<LGraph>> algorithm = lgraph.getProperty(InternalProperties.PROCESSORS);
        float monitorProgress = 1.0f / algorithm.size();
        
        if (monitor.isLoggingEnabled()) {
            // Print the algorithm configuration
            monitor.log("ELK Layered uses the following " + algorithm.size() + " modules:");
            int slot = 0;
            for (ILayoutProcessor<LGraph> processor : algorithm) {
                // SUPPRESS CHECKSTYLE NEXT MagicNumber
                String gwtDoesntSupportPrintf = (slot < 10 ? "0" : "") + (slot++);
                monitor.log("   Slot " + gwtDoesntSupportPrintf + ": " + processor.getClass().getName());
            }
        }
        
        // Invoke each layout processor
        int slotIndex = 0;
        for (ILayoutProcessor<LGraph> processor : algorithm) {
            if (monitor.isCanceled()) {
                return;
            }
            
            // Output debug graph
            // elkjs-exclude-start
            if (monitor.isLoggingEnabled()) {
                DebugUtil.logDebugGraph(monitor, lgraph, slotIndex, "Before " + processor.getClass().getSimpleName());
            }
            // elkjs-exclude-end
            
            notifyProcessorReady(lgraph, processor);
            processor.process(lgraph, monitor.subTask(monitorProgress));
            notifyProcessorFinished(lgraph, processor);
            
            slotIndex++;
        }

        // Graph debug output
        // elkjs-exclude-start
        if (monitor.isLoggingEnabled()) {
            DebugUtil.logDebugGraph(monitor, lgraph, slotIndex, "Finished");
        }
        // elkjs-exclude-end

        // Move all nodes away from the layers (we need to remove nodes from their current layer in a
        // second loop to avoid ConcurrentModificationExceptions)
        for (Layer layer : lgraph) {
            lgraph.getLayerlessNodes().addAll(layer.getNodes());
            layer.getNodes().clear();
        }
        
        for (LNode node : lgraph.getLayerlessNodes()) {
            node.setLayer(null);
        }
        
        lgraph.getLayers().clear();

        if (!monitorWasAlreadyRunning) {
            monitor.done();
        }
    }

    /**
     * Executes the given layout processor on the given list of graphs.
     *
     * @param lgraphs the list of graphs to be laid out.
     * @param monitor a progress monitor.
     * @param processor processor to execute.
     */
    private void layoutTest(final List<LGraph> lgraphs, final ILayoutProcessor<LGraph> processor) {
        // invoke the layout processor on each of the given graphs
        for (LGraph lgraph : lgraphs) {
            notifyProcessorReady(lgraph, processor);
            processor.process(lgraph, new BasicProgressMonitor());
            notifyProcessorFinished(lgraph, processor);
        }
    }


    ////////////////////////////////////////////////////////////////////////////////
    // Graph Postprocessing (Size and External Ports)

    /**
     * Sets the size of the given graph such that size constraints are adhered to.
     * Furthermore, the border spacing is added to the graph size and the graph offset.
     * Afterwards, the border spacing property is reset to 0.
     *
     * <p>Major parts of this method are adapted from
     * {@link ElkUtil#resizeNode(org.eclipse.elk.graph.ElkNode, double, double, boolean, boolean)}.</p>
     *
     * <p>Note: This method doesn't care about labels of compound nodes since those labels are not
     * attached to the graph.</p>
     *
     * @param lgraph the graph to resize.
     */
    private void resizeGraph(final LGraph lgraph) {
        Set<SizeConstraint> sizeConstraint = lgraph.getProperty(LayeredOptions.NODE_SIZE_CONSTRAINTS);
        Set<SizeOptions> sizeOptions = lgraph.getProperty(LayeredOptions.NODE_SIZE_OPTIONS);

        KVector calculatedSize = lgraph.getActualSize();
        KVector adjustedSize = new KVector(calculatedSize);

        // calculate the new size
        if (sizeConstraint.contains(SizeConstraint.MINIMUM_SIZE)) {
            KVector minSize = lgraph.getProperty(LayeredOptions.NODE_SIZE_MINIMUM);

            // if minimum width or height are not set, maybe default to default values
            if (sizeOptions.contains(SizeOptions.DEFAULT_MINIMUM_SIZE)) {
                if (minSize.x <= 0) {
                    minSize.x = ElkUtil.DEFAULT_MIN_WIDTH;
                }

                if (minSize.y <= 0) {
                    minSize.y = ElkUtil.DEFAULT_MIN_HEIGHT;
                }
            }

            // apply new size including border spacing
            adjustedSize.x = Math.max(calculatedSize.x, minSize.x);
            adjustedSize.y = Math.max(calculatedSize.y, minSize.y);
        }

        if (!lgraph.getProperty(LayeredOptions.NODE_SIZE_FIXED_GRAPH_SIZE)) {
            resizeGraphNoReallyIMeanIt(lgraph, calculatedSize, adjustedSize);
        }
    }


    /**
     * Applies a new effective size to a graph that previously had an old size calculated by the
     * layout algorithm. This method takes care of adjusting content alignments as well as external
     * ports that would be misplaced if the new size is larger than the old one.
     *
     * @param lgraph
     *            the graph to apply the size to.
     * @param oldSize
     *            old size as calculated by the layout algorithm.
     * @param newSize
     *            new size that may be larger than the old one.
     */
    private void resizeGraphNoReallyIMeanIt(final LGraph lgraph, final KVector oldSize,
            final KVector newSize) {

        // obey to specified alignment constraints
        Set<ContentAlignment> contentAlignment =
                lgraph.getProperty(LayeredOptions.CONTENT_ALIGNMENT);

        // horizontal alignment
        if (newSize.x > oldSize.x) {
            if (contentAlignment.contains(ContentAlignment.H_CENTER)) {
                lgraph.getOffset().x += (newSize.x - oldSize.x) / 2f;
            } else if (contentAlignment.contains(ContentAlignment.H_RIGHT)) {
                lgraph.getOffset().x += newSize.x - oldSize.x;
            }
        }

        // vertical alignment
        if (newSize.y > oldSize.y) {
            if (contentAlignment.contains(ContentAlignment.V_CENTER)) {
                lgraph.getOffset().y += (newSize.y - oldSize.y) / 2f;
            } else if (contentAlignment.contains(ContentAlignment.V_BOTTOM)) {
                lgraph.getOffset().y += newSize.y - oldSize.y;
            }
        }

        // correct the position of eastern and southern hierarchical ports, if necessary
        if (lgraph.getProperty(InternalProperties.GRAPH_PROPERTIES).contains(
                GraphProperties.EXTERNAL_PORTS)
                && (newSize.x > oldSize.x || newSize.y > oldSize.y)) {

            // iterate over the graph's nodes, looking for eastern / southern external ports
            // (at this point, the graph's nodes are not divided into layers anymore)
            for (LNode node : lgraph.getLayerlessNodes()) {
                // we're only looking for external port dummies
                if (node.getType() == NodeType.EXTERNAL_PORT) {
                    // check which side the external port is on
                    PortSide extPortSide = node.getProperty(InternalProperties.EXT_PORT_SIDE);
                    if (extPortSide == PortSide.EAST) {
                        node.getPosition().x += newSize.x - oldSize.x;
                    } else if (extPortSide == PortSide.SOUTH) {
                        node.getPosition().y += newSize.y - oldSize.y;
                    }
                }
            }
        }

        // Actually apply the new size
        LPadding lPadding = lgraph.getPadding();
        lgraph.getSize().x = newSize.x - lPadding.left - lPadding.right;
        lgraph.getSize().y = newSize.y - lPadding.top - lPadding.bottom;
    }

}
