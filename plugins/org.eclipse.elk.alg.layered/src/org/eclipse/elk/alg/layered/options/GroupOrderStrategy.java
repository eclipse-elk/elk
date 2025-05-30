/*******************************************************************************
 * Copyright (c) 2025 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 
 *******************************************************************************/
package org.eclipse.elk.alg.layered.options;

/**
 * Determines how to count ordering violations during layered phases.
 * Hence, this must be implemented by the several model order strategies to work correctly.
 */
public enum GroupOrderStrategy {
    /**
     * Different group are not comparable neither by their group id nor by model order. However, if a total ordering
     * is required one can of course still use either ordering to create it.
     */
    NONE,
    /**
     * The model order is more important than the group id when comparing elements from different ordering groups.
     */
    MODEL_ORDER,
    /**
     * The group id is more important than the model order then comparing elements from different ordering groups.
     */
    ENFORCED;
}
