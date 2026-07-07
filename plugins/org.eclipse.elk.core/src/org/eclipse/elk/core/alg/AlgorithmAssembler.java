/*******************************************************************************
 * Copyright (c) 2017 cds and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.alg;

import java.util.Comparator;
import java.util.EnumSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.eclipse.elk.core.util.AbstractRandomListAccessor;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

/**
 * A builder class for algorithms which are structured into phases and intermediate processing slots. The assembler
 * views layout algorithms as lists of {@link ILayoutPhase layout phases} and {@code ILayoutProcessor layout processors}
 * which, executed in order, as a pipeline, implements the algorithm.
 * 
 * <p>
 * The algorithm assembler can be used for algorithms structured into an ordered set of primary phases as defined by an
 * enumeration. Before, between, and after the phases are <em>processing slots</em> that can be filled with processors
 * that perform secondary processing on a layout graph. Processors executed in the same processing slot have to be
 * ordered in some way (see below). Different phase implementations can depend on different processors being executed,
 * which implies that the algorithm built by this assembler can completely differ subject to the chosen phase
 * implementations.
 * </p>
 * 
 * 
 * <h3>Requirements</h3>
 * 
 * <p>
 * To use this class, the following things are necessary:
 * </p>
 * <ul>
 * <li>An enumeration of the phases the algorithm is structured into. The type of this enumeration is one of the type
 * parameters of this class. This is not an enumeration of all phase implementations available for each phase.</li>
 * <li>One factory for each phase that knows how to instantiate each of the phase's available implementations. If there
 * is only a single implementation for each phase, this may be the same enumeration as the one mentioned in the previous
 * item. Usually, algorithms will offer different alternatives, which requires one enumeration for each phase.</li>
 * <li>A factory that knows how to instantiate each available layout processor. This is usually a large enumeration
 * which implements the {@link ILayoutProcessorFactory} interface.</li>
 * </ul>
 * 
 * 
 * <h3>Usage</h3>
 * 
 * <p>
 * To use this assembler, follow these steps:
 * </p>
 * <ol>
 * <li>Instantiate the algorithm assembler.</li>
 * <li>Set factories for each of the phases by calling {@link #setPhase(Enum, ILayoutPhaseFactory)}.</li>
 * <li>Optionally, call {@link #addProcessorConfiguration(LayoutProcessorConfiguration)} in order to specify layout
 * processors that should be executed in spite of no layout phase depending on them (phase-independent layout
 * processors).</li>
 * <li>Call {@link #build(Object)}.</li>
 * </ol>
 * <p>
 * A single instance can be reused by reconfiguring phases and calling {@link #build(Object)} again. Optionally, the
 * {@link #reset()} method can be called before restarting the process to reset all phase assignments and the additional
 * processing configuration.
 * </p>
 * 
 * 
 * <h3>Caching</h3>
 * 
 * <p>
 * The assembler does not use phase or processor instances themselves, but instead relies on phase and processor
 * factories. It can cache phase and processor instances obtained from these factories between subsequent calls to the
 * {@link #build(Object)} method. By default, caching is enabled.
 * </p>
 * 
 * 
 * <h3>Processor Sorting</h3>
 * 
 * <p>
 * The order of processors that share a processing slot is usually important. By default, this class uses the
 * {@link EnumBasedFactoryComparator} to sort processors by their order of appearance in their processor enumeration.
 * However, that is not the only way order can be determined. In fact, it is not even required to have all processors
 * listed in an enumeration. Order can also rely on explicitly specified dependencies. This class supports that by
 * allowing clients to specify custom comparators by calling {@link #withProcessorComparator(Comparator)}.
 * </p>
 * 
 * 
 * @param <P>
 *            enumeration of all available phases. This is not an enumeration of all phase implementations.
 * @param <G>
 *            type of the graph the created algorithm will operate on.
 * @see ILayoutProcessor
 * @see ILayoutPhase
 * @see LayoutProcessorConfiguration
 */
public final class AlgorithmAssembler<P extends Enum<P>, G>
        extends AbstractRandomListAccessor<ILayoutPhaseFactory<P, G>> {

    /** Whether layout processors and phases are cached once obtained from their factories. */
    private boolean enableCaching = true;
    /** Whether assembling the algorithm should fail if not all phases have implementations assigned. */
    private boolean failOnMissingPhase = true;
    /** The comparator used to sort processors in each processing slot. */
    private Comparator<ILayoutProcessorFactory<G>> processorComparator = new EnumBasedFactoryComparator<>();

    /** Class of the enumeration of all phases an algorithm is structured into. */
    private final Class<P> phasesEnumClass;
    /** The number of enumeration constants in the phase enumeration, saved for performance's sake. */
    private final int numberOfPhases;
    /** Set of phases whose implementations already have implementations assigned. */
    private final Set<P> configuredPhases;
    /** A configuration with additional processors which the phases don't necessarily depend upon themselves. */
    private final LayoutProcessorConfiguration<P, G> additionalProcessors;
    /** The cache for phase and processor instances, if caching is enabled. */
    private final Map<ILayoutProcessorFactory<G>, ILayoutProcessor<G>> cache;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Creation

    /**
     * Creates a new algorithm assembler based on the given enumeration of phases the algorithm should be structured
     * into.
     * 
     * @param phaseEnum
     *            enumeration class.
     * @return an algorithm assembler that can build algorithms structured into the given phases.
     */
    public static <P extends Enum<P>, G> AlgorithmAssembler<P, G> create(final Class<P> phaseEnum) {
        return new AlgorithmAssembler<>(phaseEnum);
    }

    private AlgorithmAssembler(final Class<P> phaseEnum) {
        phasesEnumClass = phaseEnum;
        numberOfPhases = phaseEnum.getEnumConstants().length;

        if (numberOfPhases == 0) {
            throw new IllegalArgumentException("There must be at least one phase in the phase enumeration.");
        }

        configuredPhases = EnumSet.noneOf(phasesEnumClass);
        additionalProcessors = LayoutProcessorConfiguration.create();
        cache = Maps.newHashMap();
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Configuration

    /**
     * Determines whether caching of layout phases and processors obtained from the factories should be enabled or not.
     * Caching is enabled by default.
     * 
     * @param enabled
     *            {@code true} if caching should be enabled.
     * @return this assembler to enable method chaining.
     */
    public AlgorithmAssembler<P, G> withCaching(final boolean enabled) {
        enableCaching = enabled;
        return this;
    }

    /**
     * Determines whether to fail upon calling the {@link #build(Object)} method if not all phases have implementations
     * configured by calling {@link #setPhase(Enum, ILayoutPhaseFactory)}. The default is to fail.
     * 
     * @param fail 
     *            {@code true} if the assembler should fail upon encountering unconfigured phases.
     * @return this assembler to enable method chaining.
     */
    public AlgorithmAssembler<P, G> withFailOnMissingPhase(final boolean fail) {
        failOnMissingPhase = fail;
        return this;
    }

    /**
     * Sets a comparator used to determine the order of processors in each processing slot. By default, processors are
     * sorted by the order of their appearance in their processor enumeration.
     * 
     * @param comparator
     *            the comparator that knows how to order processors in each slot.
     * @return this assembler to enable method chaining.
     */
    public AlgorithmAssembler<P, G> withProcessorComparator(final Comparator<ILayoutProcessorFactory<G>> comparator) {
        processorComparator = comparator;
        return this;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Phase Assembly

    /**
     * Clears the cache of layout phase and layout processor instances.
     * 
     * @return this assembler to enable method chaining.
     */
    public AlgorithmAssembler<P, G> clearCache() {
        cache.clear();
        return this;
    }

    /**
     * Resets the configured phase implementations configured through {@link #setPhase(Enum, ILayoutPhaseFactory)} and
     * any additional processors added through {@link #addProcessorConfiguration(LayoutProcessorConfiguration)}.
     * 
     * @return this assembler to enable method chaining.
     */
    public AlgorithmAssembler<P, G> reset() {
        clearList();
        configuredPhases.clear();
        additionalProcessors.clear();
        return this;
    }

    /**
     * Configures the given factory to be used to obtain an implementation for the given phase.
     * 
     * @param phase
     *            the phase to configure an implementation for.
     * @param phaseFactory
     *            the factory to use to obtain the phase implementation later.
     * @return this assembler to enable method chaining.
     */
    public AlgorithmAssembler<P, G> setPhase(final P phase, final ILayoutPhaseFactory<P, G> phaseFactory) {
        setListItem(phase.ordinal(), phaseFactory);
        configuredPhases.add(phase);
        return this;
    }

    /**
     * Adds the processors specified in the configuration to the algorithm, regardless of the requirements specified by
     * the different phases.
     * 
     * @param config
     *            a configuration that specifies additional processors to be executed.
     * @return this assembler to enable method chaining.
     */
    public AlgorithmAssembler<P, G> addProcessorConfiguration(final LayoutProcessorConfiguration<P, G> config) {
        additionalProcessors.addAll(config);
        return this;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Algorithm Building

    /**
     * Returns a list of phases and processors that, if executed in order, implement the desired algorithm. Since the
     * dependencies of phases may depend on the graph to be laid out, that graph is passed to this method.
     * 
     * @param graph
     *            the graph to be laid out.
     * @return the algorithm as a list of processors.
     * @throws IllegalStateException
     *             if {@link #failOnMissingPhase()} is active and not all phases have implementations assigned.
     */
    public List<ILayoutProcessor<G>> build(final G graph) {
        // Check if there are enough phases
        if (failOnMissingPhase && configuredPhases.size() < numberOfPhases) {
            throw new IllegalStateException("Expected " + numberOfPhases + " phases to be configured; " + "only found "
                    + configuredPhases.size());
        }

        // All the phases that could potentially have a phase implementation
        P[] phaseEnumConstants = phasesEnumClass.getEnumConstants();

        // Instantiate all configured phases
        List<ILayoutPhase<P, G>> phaseImplementations = Lists.newArrayListWithCapacity(numberOfPhases);
        for (P phase : phaseEnumConstants) {
            ILayoutPhaseFactory<P, G> phaseFactory = getListItem(phase.ordinal());
            if (phaseFactory != null) {
                phaseImplementations.add(retrievePhase(phaseFactory));
            } else {
                phaseImplementations.add(null);
            }
        }

        // Assemble a definitive processor configuration
        LayoutProcessorConfiguration<P, G> processorConfiguration = LayoutProcessorConfiguration.create();
        phaseImplementations.stream()
                .filter(phase -> phase != null)
                .map(phase -> phase.getLayoutProcessorConfiguration(graph))
                .filter(config -> config != null)
                .forEach(config -> processorConfiguration.addAll(config));
        processorConfiguration.addAll(additionalProcessors);

        // The list of processors the algorithm will be made up of
        List<ILayoutProcessor<G>> algorithm = Lists.newArrayList();

        // Add processors and phases to the algorithm
        for (P phase : phaseEnumConstants) {
            // Add processors
            algorithm.addAll(retrieveProcessors(processorConfiguration.processorsBefore(phase)));

            // Add the phase itself, if it exists
            ILayoutPhase<P, G> phaseImplementation = phaseImplementations.get(phase.ordinal());
            if (phaseImplementation != null) {
                algorithm.add(phaseImplementation);
            }
        }

        // Add processors after the last phase
        algorithm.addAll(retrieveProcessors(
                processorConfiguration.processorsAfter(phaseEnumConstants[phaseEnumConstants.length - 1])));

        return algorithm;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Utilities

    /*
     * (non-Javadoc)
     * 
     * @see org.eclipse.elk.core.alg.AbstractRandomListAccessor#provideDefault()
     */
    @Override
    protected ILayoutPhaseFactory<P, G> provideDefault() {
        return null;
    }

    /**
     * Returns a list of processors belonging to the given factories, sorted according to the configured processors
     * comparator.
     * 
     * @param factories
     *            factories for the processors to be instantiated (or retrieved from the cache).
     * @return the sorted list of instantiated processors.
     */
    private List<ILayoutProcessor<G>> retrieveProcessors(final Set<ILayoutProcessorFactory<G>> factories) {
        List<ILayoutProcessor<G>> processors = Lists.newArrayListWithCapacity(factories.size());
        factories.stream()
            .sorted(processorComparator)
            .forEach(factory -> processors.add(retrieveProcessor(factory)));
        return processors;
    }

    /**
     * Retrieves a phase obtained from the given factory either from the factory itself or from the cache. This is a
     * front for {@link #retrieveProcessor(ILayoutProcessorFactory)} that encapsulates a nasty unchecked (but safe)
     * cast.
     * 
     * @param factory
     *            the factory the phase can be obtained from.
     * @return the phase.
     */
    @SuppressWarnings("unchecked")
    private ILayoutPhase<P, G> retrievePhase(final ILayoutPhaseFactory<P, G> factory) {
        // This is an unchecked cast, but we know that the factory will return the right phase implementation since
        // the factory can only be added through the properly parameterized setPhase(...) method
        return (ILayoutPhase<P, G>) retrieveProcessor(factory);
    }

    /**
     * Returns a processor as instantiated by the given factory, either by invoking the factory or by querying the
     * cache.
     * 
     * @param factory
     *            the factory the processor can be obtained from.
     * @return the processor.
     */
    private ILayoutProcessor<G> retrieveProcessor(final ILayoutProcessorFactory<G> factory) {
        if (enableCaching) {
            if (cache.containsKey(factory)) {
                return cache.get(factory);
            } else {
                ILayoutProcessor<G> processor = factory.create();
                cache.put(factory, processor);
                return processor;
            }
        } else {
            return factory.create();
        }
    }

}
