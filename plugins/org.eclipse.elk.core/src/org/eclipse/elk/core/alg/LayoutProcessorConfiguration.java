/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
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
package org.eclipse.elk.core.alg;

import java.util.Set;

import org.eclipse.elk.core.util.AbstractRandomListAccessor;

import com.google.common.collect.Sets;

/**
 * Instances of this class specify {@link ILayoutProcessor layout processors} that should be executed in the different
 * processing slots before, in between, and after the phases of a layout algorithm. Layout processor configurations are
 * typically specified by {@link ILayoutPhase layout phases} to specify their processor dependencies. The information
 * can be used by an {@link AlgorithmAssembler} to assemble the layout algorithm. Each processor is only ever added once
 * to each processing slot, but can be part of several different slots.
 * 
 * <p>
 * Use this class by first obtaining an instance through the {@link #create()} method and then calling
 * {@link #addBefore(Enum, ILayoutProcessorFactory)} and {@link #addAfter(Enum, ILayoutProcessorFactory)} to add the
 * necessary processors before or after any of the algorithm's phases. If more than one processor needs to be added to a
 * certain processing slot, the {@link #before(Enum)} and {@link #after(Enum)} methods can be called once, followed by
 * arbitrarily many calls to {@link #add(ILayoutProcessorFactory)}.
 * </p>
 * 
 * <p>
 * Note that this class provides access to the configured processors in each slot through sets, meaning that they are
 * not sorted according to inter-processor dependencies. If used with {@link AlgorithmAssembler}, the assembler
 * determines the correct order.
 * </p>
 * 
 * @param <P>
 *            enumeration of all available phases. This is not an enumeration of all phase implementations.
 * @param <G>
 *            type of the graph the created algorithm will operate on.
 * @see ILayoutProcessor
 * @see ILayoutPhase
 * @see AlgorithmAssembler
 */
public final class LayoutProcessorConfiguration<P extends Enum<P>, G>
        extends AbstractRandomListAccessor<Set<ILayoutProcessorFactory<G>>> {

    /** The processing slot index we're currently adding processors to. */
    private int currentIndex;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Creation

    /**
     * Creates a new instance.
     * 
     * @return new layout processor configuration.
     */
    public static <P extends Enum<P>, G> LayoutProcessorConfiguration<P, G> create() {
        return new LayoutProcessorConfiguration<>();
    }
    
    /**
     * Creates a new instance which is a copy of the given instance.
     * 
     * @param source
     *            the existing configuration to copy.
     * @return new layout processor configuration.
     */
    public static <P extends Enum<P>, G> LayoutProcessorConfiguration<P, G> createFrom(
            final LayoutProcessorConfiguration<P, G> source) {
        
        return new LayoutProcessorConfiguration<P, G>().addAll(source);
    }

    private LayoutProcessorConfiguration() {
        // Le private constructor.
        clear();
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Processor Assembly

    /**
     * Resets the configuration by removing all processors.
     * 
     * @return this configuration to enable method chaining.
     */
    LayoutProcessorConfiguration<P, G> clear() {
        clearList();
        currentIndex = -1;
        return this;
    }

    /**
     * Sets things up such that subsequent calls to {@link #add(ILayoutProcessorFactory)} will add processors to the
     * processing slot right before the given phase.
     * 
     * @param phase
     *            the phase processors should be added before.
     * @return this configuration to enable method chaining.
     */
    public LayoutProcessorConfiguration<P, G> before(final P phase) {
        currentIndex = phase.ordinal();
        return this;
    }

    /**
     * Sets things up such that subsequent calls to {@link #add(ILayoutProcessorFactory)} will add processors to the
     * processing slot right after the given phase.
     * 
     * @param phase
     *            the phase processors should be added after.
     * @return this configuration to enable method chaining.
     */
    public LayoutProcessorConfiguration<P, G> after(final P phase) {
        currentIndex = phase.ordinal() + 1;
        return this;
    }

    /**
     * Adds the given processor to the current processing slot. The current processing slot is configured by calling
     * {@link #before(Enum)} or {@link #after(Enum)}. If neither of the two methods was called before or if
     * {@link #addBefore(Enum, ILayoutProcessorFactory)} or {@link #addAfter(Enum, ILayoutProcessorFactory)} was called
     * in the interim, there is no current processing slot and the call to this method fails.
     * 
     * @param processor
     *            the processor to add to the current processing slot.
     * @return this configuration to enable method chaining.
     * @throws IllegalStateException
     *             if there is no current processing slot.
     */
    public LayoutProcessorConfiguration<P, G> add(final ILayoutProcessorFactory<G> processor) {
        if (currentIndex < 0) {
            throw new IllegalStateException("Did not call before(...) or after(...) before calling add(...).");
        }

        doAdd(currentIndex, processor);
        return this;
    }

    /**
     * Adds the given processor to the processing slot right before the given phase.
     * 
     * @param phase
     *            the phase before which to add the processor.
     * @param processor
     *            the processor to add.
     * @return this configuration to enable method chaining.
     */
    public LayoutProcessorConfiguration<P, G> addBefore(final P phase, final ILayoutProcessorFactory<G> processor) {
        // Calls to add(...) after this method without calling before(...) or after(...) will usually not do
        // what people want them to, so make sure that we have no valid current index
        currentIndex = -1;

        doAdd(phase.ordinal(), processor);
        return this;
    }

    /**
     * Adds the given processor to the processing slot right after the given phase.
     * 
     * @param phase
     *            the phase after which to add the processor.
     * @param processor
     *            the processor to add.
     * @return this configuration to enable method chaining.
     */
    public LayoutProcessorConfiguration<P, G> addAfter(final P phase, final ILayoutProcessorFactory<G> processor) {
        // Calls to add(...) after this method without calling before(...) or after(...) will usually not do
        // what people want them to, so make sure that we have no valid current index
        currentIndex = -1;

        doAdd(phase.ordinal() + 1, processor);
        return this;
    }

    /**
     * Adds the given processor to the slot with the given index.
     * 
     * @param index
     *            slot index.
     * @param processor
     *            the processor to add.
     */
    private void doAdd(final int index, final ILayoutProcessorFactory<G> processor) {
        getListItem(index).add(processor);
    }

    /**
     * Adds all processors from the given configuration to this configuration.
     * 
     * @param configuration
     *            the configuration to add to this one.
     * @return this configuration to enable method chaining.
     */
    public LayoutProcessorConfiguration<P, G> addAll(final LayoutProcessorConfiguration<P, G> configuration) {
        for (int i = 0; i < configuration.getListSize(); i++) {
            getListItem(i).addAll(configuration.getListItem(i));
        }
        return this;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Processor List Building

    /**
     * Returns the set of processors configured to be executed right before the given phase.
     * 
     * @param phase
     *            the phase.
     * @return set of processors.
     */
    public Set<ILayoutProcessorFactory<G>> processorsBefore(final P phase) {
        return processors(phase.ordinal());
    }

    /**
     * Returns the set of processors configured to be executed right after the given phase.
     * 
     * @param phase
     *            the phase.
     * @return set of processors.
     */
    public Set<ILayoutProcessorFactory<G>> processorsAfter(final P phase) {
        return processors(phase.ordinal() + 1);
    }

    /**
     * Returns the set of processors in the slot with the given index.
     * 
     * @param index
     *            slot index.
     * @return set of processors.
     */
    private Set<ILayoutProcessorFactory<G>> processors(final int index) {
        return Sets.newHashSet(getListItem(index));
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Utilities

    /*
     * (non-Javadoc)
     * 
     * @see org.eclipse.elk.core.alg.AbstractRandomListAccessor#provideDefault()
     */
    @Override
    protected Set<ILayoutProcessorFactory<G>> provideDefault() {
        return Sets.newHashSet();
    }

}
