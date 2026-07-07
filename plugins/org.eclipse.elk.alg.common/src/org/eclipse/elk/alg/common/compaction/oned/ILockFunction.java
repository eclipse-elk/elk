/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.common.compaction.oned;

import org.eclipse.elk.core.options.Direction;

/**
 *
 */
@FunctionalInterface
public interface ILockFunction {
    /**
     * @param node
     * @param direction
     * @return whether the passed node may move in the passed direction.
     */
    boolean isLocked(CNode node, Direction direction);
}
